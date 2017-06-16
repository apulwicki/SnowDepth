%% Select data

load('TopoSWE.mat')
[dataSWE, dataTOPO] = ObsInCell(fullSWE.S2.input, topo_sampled);

runs            = 100;
num_excluded    = 30;

F = fieldnames(dataTOPO.G4);
for g = 1:3; glacier = options.glacier{g};
I.(glacier) = randi([1 length(dataSWE.(glacier))],[num_excluded runs]);  
end

%% Get distribution of RMSE
LRrmse = zeros(runs,3);     SKrmse = LRrmse;    RKrmse = LRrmse;

for r = 1:runs
        display(num2str(r))
for g = 1:3; glacier = options.glacier{g};
inputSWE.(glacier) = dataSWE.(glacier);
    inputSWE.(glacier)(I.(glacier)(:,r),:) = [];

inputTOPO.(glacier) = dataTOPO.(glacier);
    for i = 1:length(F)
    inputTOPO.(glacier).(F{i})(I.(glacier)(:,r),:) = [];
    end

end

% LR
    LRtest = LinearRegression(inputSWE, inputTOPO, topo_full);
    pred = SampledCell(LRtest);

    for g = 1:3; glacier = options.glacier{g};
    T.LR.(glacier) = [dataSWE.(glacier)(I.(glacier)(:,r),1), pred.(glacier)(I.(glacier)(:,r))];
    end

    for g = 1:3; glacier = options.glacier{g};
LRrmse(r,g) = sqrt(mean((T.LR.(glacier)(:,1)-T.LR.(glacier)(:,2)).^2));
    end

% % SK
%     SKtest = KrigingR_G(inputSWE);
%     pred = SampledCell(SKtest);
% 
%     for g = 1:3; glacier = options.glacier{g};
%     T.SK.(glacier) = [dataSWE.(glacier)(I,1), pred.(glacier)(I)];
%     end
% 
%     for g = 1:3; glacier = options.glacier{g};
% SKrmse(r,g) = sqrt(mean((T.SK.(glacier)(:,1)-T.SK.(glacier)(:,2)).^2));
%     end
% 
% % RK
%     RKtest = RegressionKriging(inputSWE, inputTOPO, topo_full);
%     pred = SampledCell(RKtest);
% 
%     for g = 1:3; glacier = options.glacier{g};
%     T.RK.(glacier) = [dataSWE.(glacier)(I,1), pred.(glacier)(I)];
%     end
% 
%     for g = 1:3; glacier = options.glacier{g};
% RKrmse(r,g) =     sqrt(mean((T.RK.(glacier)(:,1)-T.RK.(glacier)(:,2)).^2));
%     end

end

%% Are fullLR beta values significant?
load TopoSWE.mat 
load Full.mat

for d = 1:8
    den = options.DenOpt{d};
[ dataSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled);

for g = 1:3;        glacier = options.glacier{g};
    mu      = fullLR.(den).coeff{[8,1:7],g};

 %Get Beta dist
    X       = struct2table(TOPOdata.(glacier)); X = X{:,:}; X = [ones(length(X),1) X];
    Y       = dataSWE.(den).(glacier)(:,1);
    
    %VarCov matrix from betas
    sigmaSq = sum((Y-X*mu).^2)/(size(X,1)-size(X,2));
    VarCov  = sigmaSq*(chol(X'*X)\inv(chol(X'*X))');
    StdErr  = sqrt(diag(VarCov));     
    sig     = (abs(mu)-StdErr)>0;
    SigBeta.(den).(glacier)  = [mu, StdErr,sig];
        SigBetaDen.(glacier)(:,d) = sig;
end
end