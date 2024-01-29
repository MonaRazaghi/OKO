function [inx_kcats,inx_cmplx] = nonZero_Kcats(ecModel_batch)

start_e_mets = 1879;
start_e_rxns = size(ecModel_batch.S,2) - length(ecModel_batch.enzymes); 

pmet_inx = find(contains(ecModel_batch.metNames,'pmet'));
pmet_inx_K = pmet_inx - start_e_mets + 1;

K = ecModel_batch.S(start_e_mets:end-1,1:start_e_rxns-1);

inx = 0;
inx_c = 0;

for i = 1:size(K,2)
    inx_nnz = find(K(:,i));
    inx_k = setdiff(inx_nnz,pmet_inx_K);

    for j = 1:length(inx_k)
        inx = inx + 1;
        inx_kcats(inx,1) = inx_k(j) + start_e_mets - 1;
        inx_kcats(inx,2) = i;
    end

    if length(inx_k) > 1
        inx_c = inx_c + 1;
        inx_cmplx{inx_c,1} = inx - length(inx_k) + 1:inx;
    end
end

end