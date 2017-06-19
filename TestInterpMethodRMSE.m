%% Select data

load('TopoSWE.mat')
[dataSWE, dataTOPO] = ObsInCell(fullSWE.S2.input, topo_sampled);

runs            = 100;
num_excluded    = 200;

F = fieldnames(dataTOPO.G4);
% for g = 1:3; glacier = options.glacier{g};
% I200.(glacier) = randi([1 length(dataSWE.(glacier))],[num_excluded runs]);  
% end

%% Get distribution of RMSE
LRrmse200 = zeros(runs,3);     SKrmse200 = LRrmse;    RKrmse200 = LRrmse;

for r = 1:runs
        display(num2str(r))
for g = 1:3; glacier = options.glacier{g};
inputSWE.(glacier) = dataSWE.(glacier);
    inputSWE.(glacier)(I200.(glacier)(:,r),:) = [];

inputTOPO.(glacier) = dataTOPO.(glacier);
    for i = 1:length(F)
    inputTOPO.(glacier).(F{i})(I200.(glacier)(:,r),:) = [];
    end

end

% LR
%     LRtest = LinearRegression(inputSWE, inputTOPO, topo_full);
%     pred = SampledCell(LRtest);
% 
%     for g = 1:3; glacier = options.glacier{g};
%     T.LR.(glacier) = [dataSWE.(glacier)(I.(glacier)(:,r),1), pred.(glacier)(I.(glacier)(:,r))];
%     end
% 
%     for g = 1:3; glacier = options.glacier{g};
% LRrmse(r,g) = sqrt(mean((T.LR.(glacier)(:,1)-T.LR.(glacier)(:,2)).^2));
%     end

% SK
    SKtest = KrigingR_G(inputSWE);
    pred = SampledCell(SKtest);

    for g = 1:3; glacier = options.glacier{g};
    T200.SK.(glacier) = [dataSWE.(glacier)(I200.(glacier)(:,r),1), pred.(glacier)(I200.(glacier)(:,r))];
    end

    for g = 1:3; glacier = options.glacier{g};
SKrmse200(r,g) = sqrt(mean((T200.SK.(glacier)(:,1)-T200.SK.(glacier)(:,2)).^2));
    end

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

%% Residuals vs predicted values

load TopoSWE.mat 
load Full.mat

    yObserved  = ObsInCell(SWE, topo_sampled);

for d = 1%:8;    den = options.DenOpt{d};
    yEstimated = SampledCell(fullLR.(den));

for g = 1:3;    glacier = options.glacier{g};

    residuals.(glacier)  = yEstimated.(glacier)-yObserved(g).swe;
figure(1)    
    subplot(1,3,g)
    plot(yObserved(g).swe,residuals.(glacier),'.')
    corr(yObserved(g).swe,residuals.(glacier))^2
    title('Observed')
figure(2)    
    subplot(1,3,g)
    plot(yEstimated.(glacier),residuals.(glacier),'.') 
    title('Estimated')
    corr(yEstimated.(glacier),residuals.(glacier))^2
end
end


