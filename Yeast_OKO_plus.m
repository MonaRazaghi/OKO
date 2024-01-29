% OKO_plus

clc
clear all

% reading the model
load ecYeastGEM_batch.mat

% to avoid numerical issues
ecModel_batch.lb = 1000 * ecModel_batch.lb;
ecModel_batch.ub(find(isinf(ecModel_batch.ub))) = 1000;
ecModel_batch.ub = 1000 * ecModel_batch.ub;

num_r = size(ecModel_batch.S,2);

% finding the available kcats in the model 
[inx_Kcats,~] = nonZero_Kcats(ecModel_batch);

original_k = zeros(size(inx_Kcats,1),1);

for i = 1:size(inx_Kcats,1)
    original_k(i) = ecModel_batch.S(inx_Kcats(i,1),inx_Kcats(i,2));
end

% kcats bounds from DLKcat
T = readtable('yeast_kcat_preds.csv');
kcats_predicted = table2cell(T);

% mapping the kcats from DLKcat to the ones from the modelinx = 0;
inx = 0;
for i = 1:size(kcats_predicted,1)
    rxn = kcats_predicted{i,1};
    prot_swiss = kcats_predicted{i,2};
    if strfind(rxn,'arm')
        rxn_temp = strcat(rxn(5:end),'No');
        prot_index = find(contains(ecModel_batch.enzymes,prot_swiss));
        prot = ecModel_batch.enzGenes(prot_index);
        rxn_indices = setdiff(find(contains(ecModel_batch.rxns,rxn_temp)),find(contains(ecModel_batch.rxns,rxn)));

        inx_r = find(contains(ecModel_batch.grRules(rxn_indices),prot));
        for j = 1:length(inx_r)
            inx = inx + 1;
            inx_met = find(contains(ecModel_batch.metNames,prot_swiss));
            kcats_predicted_arm{inx,1} = ecModel_batch.rxns{rxn_indices(inx_r(j))};
            kcats_predicted_arm{inx,2} = prot_swiss;
            kcats_predicted_arm(inx,3:6) = kcats_predicted(i,3:6);

        end
    end
end

kcats_predicted(find(contains(kcats_predicted(:,1),'arm')),:) = [];

kcats_predicted = [kcats_predicted;kcats_predicted_arm];

% continue with OKO-plus

% setting the bounds in the model
inx_Kcats(:,3:5) = zeros(size(inx_Kcats,1),3);

for i = 1:size(inx_Kcats,1)
    prot = ecModel_batch.metNames{inx_Kcats(i,1)};
    prot_name = prot(6:end);
    rxn_name = ecModel_batch.rxns{inx_Kcats(i,2)};
    if find(contains(kcats_predicted(:,1),rxn_name))
        if find(contains(kcats_predicted(:,2),prot_name))
           if ~isempty(intersect(find(contains(kcats_predicted(:,1),rxn_name)),find(contains(kcats_predicted(:,2),prot_name))))
                found = intersect(find(contains(kcats_predicted(:,1),rxn_name)),find(contains(kcats_predicted(:,2),prot_name)));
                inx_Kcats(i,3) = 1;
                min_k = -1./(cell2mat(kcats_predicted(found,3)) * 3600);
                max_k = -1./(cell2mat(kcats_predicted(found,4)) * 3600);
                inx_Kcats(i,4) = min(original_k(i),min_k);
                inx_Kcats(i,5) = max(original_k(i),max_k);
            end
        end
    end
end

inx_Kcats_bounded = inx_Kcats(find(inx_Kcats(:,3) == 1),:);
inx_Kcats_bounded(:,4:5) = -1./inx_Kcats_bounded(:,4:5); 

original_k_bounded = original_k(find(inx_Kcats(:,3) == 1));

% mets to engineer
mets_to_go = xlsread('Metabolites.xlsx','Yeast','B2:B50');

% getting the biomass and production levels from OKO results
biomass_level = xlsread('Yeast_unbounded.xlsx','H1:H49');
inc_level = xlsread('Yeast_unbounded.xlsx','G1:G49');


RESULTS = cell(length(mets_to_go),13);

for met = 1:length(mets_to_go)
    
    ecModel = ecModel_batch;

    i_met = mets_to_go(met);
    
    RESULTS{met,1} = ecModel.metNames{i_met};

    % checking for available export reaction
    i_out = intersect(find(ecModel.S(i_met,:) < 0),find(contains(ecModel.rxnNames,'exchange')));

    RESULTS{met,2} = 1;
    % adding export reaction if needed
    if isempty(i_out)
        stoichCoeffList = zeros(size(ecModel.S,1),1);
        ecModel = addReaction(ecModel,'Engineering_obj','metaboliteList',{ecModel.mets{i_met}},'stoichCoeffList',[-1]);
        ecModel.lb(end) = 0;
        i_out = size(ecModel.S,2);
        RESULTS{met,2} = 0;
    end

    num_r = size(ecModel.S,2);

    % finding optimum biomass

    Vmin = ecModel.lb;
    Vmax = ecModel.ub;


    Aineq = [];
    bineq = [];

    Aeq = ecModel.S;
    beq = ecModel.b;

    max_abs_row_Aeq = max(abs(Aeq),[],2);

    for i = 1:size(Aeq,1)
        Aeq(i,:) = Aeq(i,:)/max_abs_row_Aeq(i);
        beq(i) = beq(i)/max_abs_row_Aeq(i);
    end

    f = zeros(size(Aeq,2),1);
    f(find(ecModel.c)) = -1;

    H = zeros(size(Aeq,2));

    lower = Vmin;
    upper = Vmax;

    m = struct();
    m.obj = f';
    m.Q = sparse(H);
    m.A = [sparse(Aineq); sparse(Aeq)]; % A must be sparse
    n = size(m.A, 2);
    m.vtype = repmat('C', n, 1);
    m.sense = [repmat('<',size(Aineq,1),1); repmat('=',size(Aeq,1),1)];
    m.rhs = full([bineq(:); beq(:)]); % rhs must be dense
    m.lb = lower;
    m.ub = upper;

    params = struct();
    params.FeasibilityTol=1e-6;
    params.IntFeasTol=1e-6;
    params.Method = 5;
    params.TimeLimit = 900;

    x0 = gurobi(m,params);
    opt_bio = -x0.objval;

    % finding the maximum production of the metabolite in the whole soultion
    % space
    f = zeros(size(Aeq,2),1);
    f(i_out) = -1;

    m.obj = f';

    x1 = gurobi(m,params);
    max_met = -x1.objval;

    RESULTS{met,3} = max_met;

    % finding the maximum production of the metabolite at optimum biomass
    lower(find(ecModel.c)) = 0.99 * opt_bio; % in case of infeasibility one needs to lower the fraction
    m.lb = lower;

    f = zeros(size(Aeq,2),1);
    f(i_out) = -1;

    m.obj = f';
    x2 = gurobi(m,params);

    max_met_bio = -x2.objval;

    RESULTS{met,4} = max_met_bio;

    if max_met_bio > 0

        % finding minimum enzyme usage where both biomass and the production of
        % the metabolite is at optimum

        start_e_rxns = size(ecModel.S,2) - length(ecModel.enzymes);

        lower(find(ecModel.c)) = 0.99 * opt_bio;  % in case of infeasibility one needs to lower the fraction
        lower(i_out) = 0.99 * max_met_bio;

        m.lb = lower;

        num_e = length(ecModel.enzymes);

        f = zeros(size(Aeq,2),1);
        f(start_e_rxns:start_e_rxns+num_e-1) = 1;

        m.obj = f';
        x3 = gurobi(m,params);

        e_level = x3.x(start_e_rxns:start_e_rxns+length(ecModel.enzymes)-1);

        % bounding changes in e // num_r + num_e

        beta_e = 10^(-1);

        min_e = max((1 - beta_e) * e_level,0);
        max_e = (1 + beta_e) * e_level;

        num_rxns = num_r - 1 - length(ecModel.enzymes);
        
        M = 10^6;
        
        Aineq_v_e_1 = zeros(num_e,size(Aeq,2) + num_e);
        Aineq_v_e_1(:,num_rxns + 1:num_rxns + num_e) = -eye(num_e);
        Aineq_v_e_1(:,size(Aeq,2) + 1:size(Aeq,2) + num_e) =  M * eye(num_e);
        bineq_v_e_1 = M * ones(num_e,1) - min_e;

        Aineq_v_e_2 = zeros(num_e,size(Aeq,2) + num_e);
        Aineq_v_e_2(:,num_rxns + 1:num_rxns + num_e) = eye(num_e);
        Aineq_v_e_2(:,size(Aeq,2) + 1:size(Aeq,2) + num_e) =  M * eye(num_e);
        bineq_v_e_2 = M * ones(num_e,1) + max_e;

        % v' matrix // num_r + num_e + size(inx_Kcats_bounded,1)

        A_v_prime = zeros(size(ecModel.S,1),size(inx_Kcats_bounded,1));
        for i = 1:size(inx_Kcats_bounded,1)
            A_v_prime(inx_Kcats_bounded(i,1),i) = 1;
        end

        Aeq_v = [Aeq,zeros(size(Aeq,1),num_e),A_v_prime];
        Aineq_v_e = [Aineq_v_e_1,zeros(num_e,size(inx_Kcats_bounded,1));Aineq_v_e_2,zeros(num_e,size(inx_Kcats_bounded,1))];

        bineq_v_e = [bineq_v_e_1;bineq_v_e_2];
        beq_v = beq;

        % v' and v for the magnitude of new Ks   // num_r + num_e +  2*size(inx_Kcats_bounded,1)
        start_v = num_r + num_e + 1;

        % small changes in K
        
        beta = 10^(-8);

        A1_v = zeros(size(inx_Kcats_bounded,1),size(Aeq_v,2) + size(inx_Kcats_bounded,1));
        for i = 1:size(inx_Kcats_bounded,1)
            A1_v(i,inx_Kcats_bounded(i,2)) = -beta/(1-beta);
        end

        A1_v(:,start_v:start_v+size(inx_Kcats_bounded,1)-1) = diag(1./original_k_bounded) .* eye(size(inx_Kcats_bounded,1));

        A1_v(:,start_v+size(inx_Kcats_bounded,1):start_v + 2*size(inx_Kcats_bounded,1) - 1) =  M * eye(size(inx_Kcats_bounded,1));


        A2_v = zeros(size(inx_Kcats_bounded,1),size(Aeq_v,2) + size(inx_Kcats_bounded,1));

        for i = 1:size(inx_Kcats_bounded,1)
            A2_v(i,inx_Kcats_bounded(i,2)) = -beta/(1+beta);
        end

        A2_v(:,start_v:start_v +size(inx_Kcats_bounded,1)-1) = diag(-1./original_k_bounded) .* eye(size(inx_Kcats_bounded,1));

        A2_v(:,start_v + size(inx_Kcats_bounded,1):start_v + 2*size(inx_Kcats_bounded,1) -1) =  M * eye(size(inx_Kcats_bounded,1));

        % positive 1 - delta * k >= 0

        pos = 10^(-6);
        pos_1 = 10^(-12);

        A3_v = zeros(size(inx_Kcats_bounded,1),size(Aeq_v,2) + size(inx_Kcats_bounded,1));

        for i = 1:size(inx_Kcats_bounded,1)
            A3_v(i,inx_Kcats_bounded(i,2)) = (pos_1 - 1);
        end

        A3_v(:,start_v:start_v +size(inx_Kcats_bounded,1)-1) = diag(-1./original_k_bounded) .* eye(size(inx_Kcats_bounded,1));

        % v_j = 0 -> v'_ij = 0

        MM = 10^(6);
        A4_v = zeros(size(inx_Kcats_bounded,1),size(Aeq_v,2) + size(inx_Kcats_bounded,1));

        for i = 1:size(inx_Kcats_bounded,1)
            A4_v(i,inx_Kcats_bounded(i,2)) = -MM;
        end

        A4_v(:,start_v:start_v +size(inx_Kcats_bounded,1)-1) = -eye(size(inx_Kcats_bounded,1));


        A5_v = zeros(size(inx_Kcats_bounded,1),size(Aeq_v,2) + size(inx_Kcats_bounded,1));

        for i = 1:size(inx_Kcats_bounded,1)
            A5_v(i,inx_Kcats_bounded(i,2)) = -MM;
        end

        A5_v(:,start_v:start_v +size(inx_Kcats_bounded,1)-1) = eye(size(inx_Kcats_bounded,1));

        % controling the magnitude of changed Kcats
        
        A6_v = zeros(size(inx_Kcats_bounded,1),size(Aeq_v,2) + size(inx_Kcats_bounded,1));

        for i = 1:size(inx_Kcats_bounded,1)
            A6_v(i,inx_Kcats_bounded(i,2)) = inx_Kcats_bounded(i,4) + 1/original_k_bounded(i);
        end

        A6_v(:,start_v:start_v +size(inx_Kcats_bounded,1)-1) = diag((1./original_k_bounded).*inx_Kcats_bounded(:,4));

        %
        A7_v = zeros(size(inx_Kcats_bounded,1),size(Aeq_v,2) + size(inx_Kcats_bounded,1));

        for i = 1:size(inx_Kcats_bounded,1)
            A7_v(i,inx_Kcats_bounded(i,2)) = -inx_Kcats_bounded(i,5) - 1/original_k_bounded(i);
        end

        A7_v(:,start_v:start_v +size(inx_Kcats_bounded,1)-1) = diag((-1./original_k_bounded).*inx_Kcats_bounded(:,5));

        %
        Aineq_final = [Aineq_v_e,zeros(size(Aineq_v_e,1),size(inx_Kcats_bounded,1));A1_v;A2_v;A3_v;A4_v;A5_v;A6_v;A7_v];
        bineq_final = [bineq_v_e;M*ones(2*size(inx_Kcats_bounded,1),1);zeros(size(inx_Kcats_bounded,1),1);pos*ones(2*size(inx_Kcats_bounded,1),1);zeros(2*size(inx_Kcats_bounded,1),1)];

        Aeq_final = [Aeq_v,zeros(size(Aeq_v,1),size(inx_Kcats_bounded,1))];
        beq_final = beq_v;

       
        lower = [lower;zeros(num_e,1);-10^(12)*ones(size(inx_Kcats_bounded,1),1);zeros(size(inx_Kcats_bounded,1),1)];
        upper = [upper;ones(num_e,1);10^(12)*ones(size(inx_Kcats_bounded,1),1);ones(size(inx_Kcats_bounded,1),1)];

        int = union(num_r + 1:num_r+num_e,start_v+size(inx_Kcats_bounded,1):start_v + 2*size(inx_Kcats_bounded,1) - 1);

        %% Engineering

        lower(i_out) = 0;
        % product level from OKO results
        upper(i_out) = (inc_level(met)/0.99) * max_met_bio;
        % biomass level from OKO results
        lower(find(ecModel.c)) =  (biomass_level(met))/0.99;

        f = zeros(num_r + num_e + 2 * size(inx_Kcats_bounded,1),1);
        f(i_out) = -1;

        H = zeros(size(Aeq_final,2));

        m_eng_temp = struct();
        m_eng_temp.obj = f';
        m_eng_temp.Q = sparse(H);
        m_eng_temp.A = [sparse(Aineq_final); sparse(Aeq_final)]; % A must be sparse
        n = size(m_eng_temp.A, 2);
        m_eng_temp.vtype = repmat('C', n, 1);
        m_eng_temp.vtype(int) = 'I';
        m_eng_temp.sense = [repmat('<',size(Aineq_final,1),1); repmat('=',size(Aeq_final,1),1)];
        m_eng_temp.rhs = full([bineq_final(:); beq_final(:)]); % rhs must be dense
        m_eng_temp.lb = lower;
        m_eng_temp.ub = upper;

        x_eng_temp = gurobi(m_eng_temp,params);

        % in case of infeasibilty we remove the upper bound for product

        if strcmp(x_eng_temp.status,'INFEASIBLE')
            warning('Second objective')
            lower(i_out) = (inc_level(met)/0.99) * max_met_bio;
            upper(i_out) = 10^6;
            lower(find(ecModel.c)) =  (biomass_level(met))/0.99;

            m_eng_temp.lb = lower;
            m_eng_temp.ub = upper;

            f = zeros(num_r + num_e + 2 * size(inx_Kcats_bounded,1),1);
            f(i_out) = -1;

            m_eng_temp.obj = f';

            x_eng_temp = gurobi(m_eng_temp,params);
        end

        RESULTS{met,5} = x_eng_temp.x(i_out);

        lower(i_out) = 0.99 * x_eng_temp.x(i_out);
        lower(find(ecModel.c)) = 0.99 * x_eng_temp.x(find(ecModel.c));
         
        sol = x_eng_temp.x;
        sol(find(sol(1:num_r) < 0)) = 0;
        i_lower = find(sol < lower);
        sol(i_lower) = lower(i_lower);

        i_upper = find(sol > upper);
        sol(i_upper) = upper(i_upper);

        f = zeros(num_r + num_e + 2 * size(inx_Kcats_bounded,1),1);
        f(num_r + 1:num_r + num_e) = -10;
        f(num_r + num_e + size(inx_Kcats_bounded,1) + 1:num_r + num_e + 2 * size(inx_Kcats_bounded,1)) = -1;

        H = zeros(size(Aeq_final,2));

        m_eng = struct();
        m_eng.obj = f';
        m_eng.Q = sparse(H);
        m_eng.A = [sparse(Aineq_final); sparse(Aeq_final)]; % A must be sparse
        n = size(m_eng.A, 2);
        m_eng.vtype = repmat('C', n, 1);
        m_eng.vtype(int) = 'I';
        m_eng.sense = [repmat('<',size(Aineq_final,1),1); repmat('=',size(Aeq_final,1),1)];
        m_eng.rhs = full([bineq_final(:); beq_final(:)]); % rhs must be dense
        m_eng.lb = lower;
        m_eng.ub = upper;
        m_eng.start = sol;

        x_eng = gurobi(m_eng,params);

        RESULTS{met,6} = x_eng.x(i_out);
        RESULTS{met,7} = x_eng.x(i_out)/max_met_bio;
        RESULTS{met,8} = x_eng.x(find(ecModel.c));

        v_prime = x_eng.x(num_r + num_e + 1:num_r + num_e + size(inx_Kcats_bounded,1));
        y_k = x_eng.x(num_r + num_e + size(inx_Kcats_bounded,1) + 1:num_r + num_e + 2*size(inx_Kcats_bounded,1));
        y_e = x_eng.x(num_r + 1:num_r + num_e);

        delta_e = zeros(size(inx_Kcats_bounded,1),1);
        v = zeros(size(inx_Kcats_bounded,1),1);
        k_new = zeros(size(inx_Kcats_bounded,1),1);
        
        count = 0; 
        for i = 1:size(inx_Kcats_bounded,1)
            inx_j = inx_Kcats_bounded(i,2);
            v(i) = x_eng.x(inx_j);

            delta_e(i) = v_prime(i)/v(i);
            k = (-1/original_k_bounded(i));
            k_new(i) = k/(1-delta_e(i)*k);
        end


        inx_changed_delta = find(y_k == 0);
        inx_changed_e = find(y_e == 0) + num_rxns;

        inx_changed_delta = setdiff(inx_changed_delta,union(find(delta_e == 0),union(find(isinf(delta_e)),find(isnan(delta_e)))));

        sum_change_k = 0;
        dec_k = 0;
        inc_k = 0;

        dec_e = 0;
        inc_e = 0;

        Results_delta = cell(length(inx_changed_delta),7);

        for i = 1:length(inx_changed_delta)
            inx = inx_changed_delta(i);

            if (1 - delta_e(inx)*(-1/original_k(inx))) < 10^(-6)
                continue
            end

            enz = ecModel.metNames{inx_Kcats_bounded(inx,1)};
            Results_delta{i,1} = enz(6:end);
            Results_delta{i,2} = ecModel.rxnNames{inx_Kcats_bounded(inx,2)};
            Results_delta{i,3} = -1/original_k_bounded(inx);
            Results_delta{i,4} = Results_delta{i,3}/(1 - delta_e(inx)*Results_delta{i,3});
            Results_delta{i,5} = inx_Kcats_bounded(inx,4);
            Results_delta{i,6} = inx_Kcats_bounded(inx,5);
            Results_delta{i,7} = log(Results_delta{i,4}/Results_delta{i,3});
            
            if Results_delta{i,5} > 0
                inc_k = inc_k + 1;
            end

            if Results_delta{i,5} < 0
                dec_k = dec_k + 1;
            end

            sum_change_k = sum_change_k + abs((delta_e(inx) * Results_delta{i,3} * Results_delta{i,3})/(1 - delta_e(inx)*Results_delta{i,3}))/Results_delta{i,3};
        end

        for i = 1:length(inx_changed_e)
            inx = inx_changed_e(i);

            if x_eng.x(inx) < 10^(-6)
                if max_e(inx - num_rxns) < 10^(-6)
                    continue
                end
            end


            enz = ecModel.rxnNames{inx};
            Results_delta{i + length(inx_changed_delta),1} = enz(11:end);
            Results_delta{i+ length(inx_changed_delta),2} = x_eng.x(inx);
            Results_delta{i+ length(inx_changed_delta),3} = min_e(inx - num_rxns);
            Results_delta{i+ length(inx_changed_delta),4} = max_e(inx - num_rxns);
                        
            if Results_delta{i+ length(inx_changed_delta),2} < Results_delta{i+ length(inx_changed_delta),3}
                dec_e = dec_e + 1;
            end

            if Results_delta{i+ length(inx_changed_delta),2} > Results_delta{i+ length(inx_changed_delta),4}
                inc_e = inc_e + 1;
            end

        end

        RESULTS{met,9} = inc_k;
        RESULTS{met,10} = dec_k;
        RESULTS{met,11} = sum_change_k;

        RESULTS{met,12} = inc_e;
        RESULTS{met,13} = dec_e;

        row_has_empty = all(cellfun(@isempty, Results_delta), 2);
        Results_delta(row_has_empty,:) = [];

        name = strcat(strcat('met_',num2str(met)),'_OKO_plus.txt');
        writecell(Results_delta,name,'Delimiter','tab');

        clear Aineq_final
        clear bineq_final
        clear lower
        clear upper
        clear int
        clear x0
        clear x1
        clear x2
        clear x3
        clear x_eng_temp
        clear x_eng_temp_2
        clear x_eng

    end
end
