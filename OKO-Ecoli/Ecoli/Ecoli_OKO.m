% OKO

clc
clear all

% reading the model
load eciML1515_batch.mat

% to avoid numerical issues
ecModel_batch.lb = 1000 * ecModel_batch.lb;
ecModel_batch.ub = 1000 * ecModel_batch.ub;

num_r = size(ecModel_batch.S,2);

% finding the available kcats in the model 
[inx_Kcats,inx_cmplx] = nonZero_Kcats(ecModel_batch);


original_k = zeros(size(inx_Kcats,1),1);

for i = 1:size(inx_Kcats,1)
    original_k(i) = ecModel_batch.S(inx_Kcats(i,1),inx_Kcats(i,2));
end


% mets to engineer
mets_to_go = xlsread('Metabolites.xlsx','Ecoli','B2:B42');

RESULTS = cell(length(mets_to_go),10);

for met = 1:length(mets_to_go)

    ecModel = ecModel_batch;

    i_met = mets_to_go(met);
    RESULTS{met,1} = ecModel.metNames{i_met};

    % checking for available export reaction
    i_out = intersect(find(ecModel.S(i_met,:) < 0),find(contains(ecModel.rxns,'EX')));
        
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

    lower = ecModel.lb;
    upper = ecModel.ub;

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

        lower(find(ecModel.c)) = 0.98 * opt_bio; % in case of infeasibility one needs to lower the fraction
        lower(i_out) = 0.99 * max_met_bio;

        m.lb = lower;

        num_e = length(ecModel.enzymes);

        f = zeros(size(Aeq,2),1);
        f(start_e_rxns:start_e_rxns+num_e-1) = 1;

        m.obj = f';
        x3 = gurobi(m,params);

        e_level = x3.x(start_e_rxns:start_e_rxns+length(ecModel.enzymes)-1);

        % bounding changes in e // num_r 

        beta_e = 10^(-1);

        min_e = max((1 - beta_e) * e_level,0);
        max_e = (1 + beta_e) * e_level;

        num_rxns = num_r - 1 - length(ecModel.enzymes);

        lower(num_rxns + 1:num_rxns + num_e) = min_e;
        upper(num_rxns + 1:num_rxns + num_e) = max_e;


        % v' matrix // num_r + size(inx_Kcats,1)


        A_v_prime = zeros(size(ecModel.S,1),size(inx_Kcats,1));
        for i = 1:size(inx_Kcats,1)
            A_v_prime(inx_Kcats(i,1),i) = 1;
        end

        Aeq_v = [Aeq,A_v_prime];
        Aineq_v = [];

        bineq_v = bineq;
        beq_v = beq;

        % v' for complexes
        start_v = num_r + 1;

        for i = 1:length(inx_cmplx)
            inx_v_primes = inx_cmplx{i};
            for j = 1:length(inx_v_primes)-1
                Aeq_v(end+1,inx_v_primes(j) + start_v - 1) = 1;
                Aeq_v(end,inx_v_primes(j+1) + start_v - 1) = -1;
                beq_v(end+1) = 0;
            end
        end

        % v' and v for the magnitude of new Ks   // num_r + 2*size(inx_Kcats,1)
        % small changes in K
        M = 10^6;
        beta = 10^(-8);

        A1_v = zeros(size(inx_Kcats,1),size(Aeq_v,2) + size(inx_Kcats,1));
        for i = 1:size(inx_Kcats,1)
            A1_v(i,inx_Kcats(i,2)) = -beta/(1-beta);
        end

        A1_v(:,start_v:start_v+size(inx_Kcats,1)-1) = diag(1./original_k) .* eye(size(inx_Kcats,1));

        A1_v(:,start_v+size(inx_Kcats,1):start_v + 2*size(inx_Kcats,1) - 1) =  M * eye(size(inx_Kcats,1));


        A2_v = zeros(size(inx_Kcats,1),size(Aeq_v,2) + size(inx_Kcats,1));

        for i = 1:size(inx_Kcats,1)
            A2_v(i,inx_Kcats(i,2)) = -beta/(1+beta);
        end

        A2_v(:,start_v:start_v +size(inx_Kcats,1)-1) = diag(-1./original_k) .* eye(size(inx_Kcats,1));

        A2_v(:,start_v + size(inx_Kcats,1):start_v + 2*size(inx_Kcats,1) -1) =  M * eye(size(inx_Kcats,1));

        % positive 1 - delta * k >= 0
        pos = 10^(-6);
        pos_1 = 10^(-12);

        A3_v = zeros(size(inx_Kcats,1),size(Aeq_v,2) + size(inx_Kcats,1));

        for i = 1:size(inx_Kcats,1)
            A3_v(i,inx_Kcats(i,2)) = (pos_1 - 1);
        end

        A3_v(:,start_v:start_v +size(inx_Kcats,1)-1) = diag(-1./original_k) .* eye(size(inx_Kcats,1));

        % v_j = 0 -> v'_ij = 0

        MM = 10^(6);
        A4_v = zeros(size(inx_Kcats,1),size(Aeq_v,2) + size(inx_Kcats,1));

        for i = 1:size(inx_Kcats,1)
            A4_v(i,inx_Kcats(i,2)) = -MM;
        end

        A4_v(:,start_v:start_v +size(inx_Kcats,1)-1) = -eye(size(inx_Kcats,1));


        A5_v = zeros(size(inx_Kcats,1),size(Aeq_v,2) + size(inx_Kcats,1));

        for i = 1:size(inx_Kcats,1)
            A5_v(i,inx_Kcats(i,2)) = -MM;
        end

        A5_v(:,start_v:start_v +size(inx_Kcats,1)-1) = eye(size(inx_Kcats,1));

        % controling the magnitude of changed Kcats
        alpha = 10^(1);

        A6_v = zeros(size(inx_Kcats,1),size(Aeq_v,2) + size(inx_Kcats,1));

        for i = 1:size(inx_Kcats,1)
            A6_v(i,inx_Kcats(i,2)) = (1 - alpha);
        end

        A6_v(:,start_v:start_v +size(inx_Kcats,1)-1) = diag(1./original_k);

        %
        A7_v = zeros(size(inx_Kcats,1),size(Aeq_v,2) + size(inx_Kcats,1));

        for i = 1:size(inx_Kcats,1)
            A7_v(i,inx_Kcats(i,2)) = (1/alpha - 1);
        end

        A7_v(:,start_v:start_v +size(inx_Kcats,1)-1) = diag(-1./original_k);

        %
        Aineq_final = [A1_v;A2_v;A3_v;A4_v;A5_v;A6_v;A7_v];
        bineq_final = [M*ones(2*size(inx_Kcats,1),1);zeros(size(inx_Kcats,1),1);pos*ones(2*size(inx_Kcats,1),1);zeros(2*size(inx_Kcats,1),1)];

        Aeq_final = [Aeq_v,zeros(size(Aeq_v,1),size(inx_Kcats,1))];
        beq_final = beq_v;


        lower = [lower;-10^(12)*ones(size(inx_Kcats,1),1);zeros(size(inx_Kcats,1),1)];
        upper = [upper;10^(12)*ones(size(inx_Kcats,1),1);ones(size(inx_Kcats,1),1)];

        int = start_v+size(inx_Kcats,1):start_v + 2*size(inx_Kcats,1) - 1;



        %% Engineering

        % initial solution
        fake_sol = zeros(size(Aineq_final,2),1);
        x3.x(find(x2.x(1:num_r) < 0)) = 0;
        fake_sol(1:num_r) = x3.x;
        fake_sol(start_v:start_v + size(inx_Kcats,1) - 1) = 0;
        fake_sol(start_v + size(inx_Kcats,1):start_v + 2*size(inx_Kcats,1) - 1) = 1;
        
        % we first put an upper bound for the biomass to make it almost the
        % same as what we had in the wild-type
        lower(i_out) = 0;
        upper(i_out) = 1.98 * max_met_bio;
        lower(find(ecModel.c)) = 0.5 * opt_bio;
        upper(find(ecModel.c)) =  0.99 * opt_bio;

        f = zeros(num_r + 2 * size(inx_Kcats,1),1);
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
        m_eng_temp.start = fake_sol;

        x_eng_temp = gurobi(m_eng_temp,params);

        % in case of infeasibilty we remove the upper bound for biomass
        if strcmp(x_eng_temp.status,'INFEASIBLE')
            lower(find(ecModel.c)) =  0.5 * opt_bio;
            upper(find(ecModel.c)) = 10^6;

            m_eng_temp.lb = lower;
            m_eng_temp.ub = upper;

            f = zeros(num_r + 2 * size(inx_Kcats,1),1);
            f(i_out) = -1;

            m_eng_temp.obj = f';

            x_eng_temp = gurobi(m_eng_temp,params);
        end

        RESULTS{met,5} = x_eng_temp.x(i_out);
        lower(i_out) = 0.9 * x_eng_temp.x(i_out); % in case of infeasibility one needs to lower the fraction
        lower(find(ecModel.c)) = 0.9 * x_eng_temp.x(find(ecModel.c)); % in case of infeasibility one needs to lower the fraction

        sol = x_eng_temp.x;
        sol(find(sol(1:num_r) < 0)) = 0;
        i_lower = find(sol < lower); 
        sol(i_lower) = lower(i_lower);

        i_upper = find(sol > upper);
        sol(i_upper) = upper(i_upper);

        f = zeros(num_r + 2 * size(inx_Kcats,1),1);
        f(num_r + size(inx_Kcats,1) + 1:num_r + 2 * size(inx_Kcats,1)) = -1;

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

        RESULTS{met,6} = x_eng.x(i_out); % production level in the engineered model
        RESULTS{met,7} = RESULTS{met,6}/max_met_bio; % fold change of the production
        RESULTS{met,8} = x_eng.x(find(ecModel.c)); % biomass level in the engineered model

        % finding the modifications
        v_prime = x_eng.x(num_r + 1:num_r + size(inx_Kcats,1));
        y_k = x_eng.x(num_r + size(inx_Kcats,1) + 1:num_r + 2*size(inx_Kcats,1));

        delta_e = zeros(size(inx_Kcats,1),1);
        v = zeros(size(inx_Kcats,1),1);

        for i = 1:size(inx_Kcats,1)
            inx_j = inx_Kcats(i,2);
            v(i) = x_eng.x(inx_j);

            delta_e(i) = v_prime(i)/v(i);
        end


        inx_changed_delta = find(y_k == 0);

        inx_changed_delta = setdiff(inx_changed_delta,union(find(delta_e == 0),union(find(isinf(delta_e)),find(isnan(delta_e)))));

        sum_change_k = 0;
        dec_k = 0;
        inc_k = 0;

        Results_delta = cell(length(inx_changed_delta),6);

        for i = 1:length(inx_changed_delta)
            inx = inx_changed_delta(i);
    
            if (1 - delta_e(inx)*(-1/original_k(inx))) < 10^(-6)
                continue
            end

            enz = ecModel.metNames{inx_Kcats(inx,1)};
            Results_delta{i,1} = enz(6:end);
            Results_delta{i,2} = ecModel.rxnNames{inx_Kcats(inx,2)};
            Results_delta{i,3} = -1/original_k(inx);
            Results_delta{i,4} = Results_delta{i,3}/(1 - delta_e(inx)*Results_delta{i,3});
            Results_delta{i,5} = (delta_e(inx) * Results_delta{i,3} * Results_delta{i,3})/(1 - delta_e(inx)*Results_delta{i,3});
            Results_delta{i,6} = log(Results_delta{i,4}/Results_delta{i,3});

            if Results_delta{i,5} > 0
                inc_k = inc_k + 1;
            end

            if Results_delta{i,5} < 0
                dec_k = dec_k + 1;
            end

            sum_change_k = sum_change_k + abs((delta_e(inx) * Results_delta{i,3} * Results_delta{i,3})/(1 - delta_e(inx)*Results_delta{i,3}))/Results_delta{i,3};
        end
        
        RESULTS{met,9} = inc_k; %number of increased kcats
        RESULTS{met,10} = dec_k; %number of decreased kcats
        RESULTS{met,11} = sum_change_k;

        row_has_empty = any(cellfun(@isempty, Results_delta), 2);
        Results_delta(row_has_empty,:) = [];

        name = strcat(strcat('met_',num2str(met)),'_OKO.txt');
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
