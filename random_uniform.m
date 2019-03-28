function [swe_out, topo_out] = random_uniform(subsetSWE, subsetTOPO, n)

    E = subsetSWE(:,2);
    N = subsetSWE(:,3);
    SWE = subsetSWE(:,1);
    
    m = length(N);
    dmetric = 0;
    trials = 5;

for k=1:trials
    pp=randperm(m);
    Nset = N(pp(1:n));
    Eset = E(pp(1:n));
    SWEset = SWE(pp(1:n));
    
%     figure(1)
%     plot(Eset, Nset, '.')
%     pause
    
    dset = zeros(n,n);
    for i=1:n
        for j=1:n
            dset(i,j) = sqrt((Nset(i)-Nset(j)).^2 + (Eset(i) - Eset(j)).^2);
        end
    end
    
     dbest = min(min(dset + 1e6*eye(n,n)));
    if dbest > dmetric
%            display([k, dbest, dmetric, SWEset(1)]);
        dmetric = dbest;
        Ebest = Eset;
        Nbest = Nset;
        SWEbest = SWEset;
        
        param = fieldnames(subsetTOPO);
        for t = 1:length(param)
           TOPObest.(param{t}) = subsetTOPO.(param{t})(pp(1:n));
        end
    end
end
    swe_out = [SWEbest, Ebest, Nbest];
    topo_out = TOPObest;
end
