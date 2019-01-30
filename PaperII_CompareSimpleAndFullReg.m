%% Compare beta coefficients and resulting distributions from simple LR and BMA

load PaperII_Syntheic.mat
load PaperII_Compare_allbeta.mat snowdist_model coeffSYN

file_path = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/';
%file_path = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/';

noise_number = 100;
num_models = 10;
nn = 18;
%% Get synthetic snow distributions

% Beta coeff ranges
elev_range  = [0.585, 0.810];
dc_range    = [0, 0.015];
aspect_range = [0, 0.093];
slope_range = [-0.09, 0.277];
N_range     = [-0.414, 0.172];
curv_range  = [-0.077, 0.02];
sx_range    = [-0.294, 0.260];

intercept   = [0.6205, 0.2627, 0.2354];

% Normal random beta values

elev_beta   = normrnd(mean(elev_range), std(elev_range), [1,num_models]);
dc_beta     = normrnd(mean(dc_range), std(dc_range), [1,num_models]);
aspect_beta = normrnd(mean(aspect_range), std(aspect_range), [1,num_models]);
slope_beta  = normrnd(mean(slope_range), std(slope_range), [1,num_models]);
N_beta      = normrnd(mean(N_range),std(N_range),[1,num_models]);
curv_beta   = normrnd(mean(curv_range), std(curv_range), [1,num_models]);
sx_beta     = normrnd(mean(sx_range),   std(sx_range),   [1,num_models]);

coeffSYN = [elev_beta', dc_beta', aspect_beta', slope_beta', N_beta', curv_beta', sx_beta'];
% Generate snow dist models
% load TopoSWE.mat topo_full
% run OPTIONS.m

    % Remove dc, aspect and Northness
%     for g = 1:3;    glacier = options.glacier{g};
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
%     end

for m = 1:num_models

for g = 1:3;    glacier = options.glacier{g};
    syn_model   = repmat(intercept(g), options.mapsize(g,:));
%     betaCoeff   = [elev_beta(m), sx_beta(m), curv_beta(m), slope_beta(m)];   
    betaCoeff   = coeffSYN(m,:);   
    topoCoeff   = fieldnames(topo_full.G4);
    
    for c = 1:length(betaCoeff)
        param     = topoCoeff{c};
        sweT      = topo_full.(glacier).(param)*betaCoeff(c);
        syn_model = syn_model + sweT;
    end
        
    syn_model(syn_model<0) = 0; %Set min to 0
   
    snowdist_model(m).(glacier) = syn_model;
    
end
end
% save('PaperII_Compare_allbeta.mat','snowdist_model','coeffSYN')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Plot synthetic models
%     load PaperII_Syntheic.mat
%     load Full.mat fullSWE
% 
% for m =1:num_models
%     figure(m);
%     PlotTopoParameter(snowdist_model(m), 'none', 'WB (m w.e.)', fullSWE.S2.input, 'black', 'none')
% end

%% SAMPLING THEORETICAL FIELD FROM PATTERNS

% Get WB field that is the "true" field
   % clear; close all
% load TopoSWE.mat topo_full

 %Make matrix with utm of each grid cell
for g = 1:3;    glacier = options.glacier{g};
minE = min(options.rig.(glacier)(:,1));
minN = min(options.rig.(glacier)(:,2));

    nE = options.mapsize(g,2);            
    nN = options.mapsize(g,1);
    
    utmGridE.(glacier) = repmat([1:nE]*40+minE,nN,1);   
    utmGridN.(glacier) = repmat([nN:-1:1]'*40+minN,1,nE); 
end
%  clear g* min* n*

for mc = 1:num_models

%Additing noise to distributed WB
fullWB = snowdist_model(mc);
 
% Obtaining pattern data for utm, topo, and wb

 %Get csv files
pattern.Circle = csvread([file_path,'CellNum_Circle.csv'],1,0);
pattern.Centreline = csvread([file_path,'CellNum_Centreline.csv'],1,0);
pattern.CentreTransect = csvread([file_path,'CellNum_Transverse.csv'],1,0);
pattern.Hourglass = csvread([file_path,'CellNum_Hourglass.csv'],1,0);
pattern.HourCircle = csvread([file_path,'CellNum_HourglassCircle.csv'],1,0);

namesP = fieldnames(pattern);
for p = 1:length(namesP)
        [~,ia,~] = unique(pattern.(namesP{p})(:,3));
        ia = sort(ia);
    pattern.(namesP{p}) = pattern.(namesP{p})(ia,:);
    
 clear I*
for g = 1:3;    glacier = options.glacier{g};
        I_E(:,1) = pattern.(namesP{p})(:,1) > min(utmGridE.(glacier)(:));
        I_E(:,2) = pattern.(namesP{p})(:,1) < max(utmGridE.(glacier)(:));
    I_E = all(I_E,2);
        I_N(:,1) = pattern.(namesP{p})(:,2) > min(utmGridN.(glacier)(:));
        I_N(:,2) = pattern.(namesP{p})(:,2) < max(utmGridN.(glacier)(:));
    I_N = all(I_N,2);
    I = all([I_E I_N],2);
    
    utmPattern.(namesP{p}).(glacier)(:,1:3) = pattern.(namesP{p})(I,:); 
end

 %Get UTM and WB for cells at the pattern locations
 clear e n wb
 for g = 1:3;    glacier = options.glacier{g};
    e = zeros(length(utmPattern.(namesP{p}).(glacier)),1); num_models = e; wb = e;

    for i = 1:length(utmPattern.(namesP{p}).(glacier))
        [~, t1] =  min(abs(utmGridE.(glacier)(1,:) - utmPattern.(namesP{p}).(glacier)(i,1)));
       pUTM.(namesP{p}).(glacier)(i,1) = utmGridE.(glacier)(1,t1);
            [~, t2] =  min(abs(utmGridN.(glacier)(:,1) - utmPattern.(namesP{p}).(glacier)(i,2))); 
       pUTM.(namesP{p}).(glacier)(i,2) = utmGridN.(glacier)(t2,1); 
       
       pWB.(namesP{p}).(glacier)(i,1) = fullWB.(glacier)(t2,t1);
       
            topoparams = fieldnames(topo_full.G4);
       for f = 1:length(topoparams)
       pTOPO.(namesP{p}).(glacier).(topoparams{f})(i,1) = topo_full.(glacier).(topoparams{f})(t2,t1); 
       end
    end
end
end

% RANDOM - SAFE AREA   
sizeR = 200;

cells = dlmread([file_path,'SafeArea.csv']);

cells = cells(1:5:end-5,2:3);   
Gsplit = [1, 1950, 3235, length(cells)]; Gsplit = flip(Gsplit);
    for g = 1:3;    glacier = options.glacier{g};
        safeR.(glacier) = cells(Gsplit(g+1):Gsplit(g)-1,:);
            I = randi(length(safeR.(glacier)), sizeR, 1); 
        safeR.(glacier) = safeR.(glacier)(I,:);

    end

 %Get UTM and WB for cells at the pattern locations
 clear e n wb
 for g = 1:3;    glacier = options.glacier{g};
    e = zeros(length(safeR.(glacier)),1); num_models = e; wb = e;

    for i = 1:length(safeR.(glacier))
        [~, t1] =  min(abs(utmGridE.(glacier)(1,:) - safeR.(glacier)(i,1)));
       pUTM.RandomSafe.(glacier)(i,1) = utmGridE.(glacier)(1,t1);
            [~, t2] =  min(abs(utmGridN.(glacier)(:,1) - safeR.(glacier)(i,2))); 
       pUTM.RandomSafe.(glacier)(i,2) = utmGridN.(glacier)(t2,1); 
       
       pWB.RandomSafe.(glacier)(i,1) = fullWB.(glacier)(t2,t1);
       
            topoparams = fieldnames(topo_full.G4);
       for f = 1:length(topoparams)
       pTOPO.RandomSafe.(glacier).(topoparams{f})(i,1) = topo_full.(glacier).(topoparams{f})(t2,t1); 
       end
    end
 end

% load PaperII_AblationArea.mat AblationArea
for g=1:3; glacier = options.glacier{g};
AblationArea.(glacier)(AblationArea.(glacier)==-0.1)=NaN;
AblationArea.(glacier)(~isnan(AblationArea.(glacier)))=1;
end
    
%%
% namesP = fieldnames(pWB);
 namesP = {'CentreTransect'};

real_measure = SampledCell(snowdist_model(mc));

for p = 1:length(namesP)

for ss = nn %6:45
   [WBinput, TOPOinput, UTMinput] = SubsetSampleSize( pWB, pTOPO, pUTM, ss );

       display([' Sample size: ',num2str(ss),' Pattern: ',namesP{p}, ' Run:',num2str(mc)])
    clear syn_measure_tmp syn_measure real_measure_tmp
    %Add some noise
    WBinputN = WBnoise(WBinput.(namesP{p}),'low');
    
        % BMA
        swe_input.(glacier) = WBinputN.(glacier);
        topo_input = TOPOinput.(namesP{p});

        cd BMS
        [BMSinit, BMSres] = BMS_R(swe_input, topo_input);
        cd ..

        for g = 3;        glacier = char(options.glacier(g));
        BMS.(glacier) = BMSinit.(glacier)(:,1);   
        coeffsBMA(ss).(namesP{p})(mc).(glacier)   = BMS.(glacier){1:8,1};

        swe_pred = repmat(BMS.(glacier){8,1}, options.mapsize(g,:));
        betaCoeff = BMS.(glacier){1:7,1};    topoCoeff = fieldnames(topo_full.G4);
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;
        swe_pred = swe_pred.*AblationArea.(glacier); 
        predBMA(ss).(namesP{p})(mc).(glacier) = swe_pred;
        
        sampledtemp             = swe_pred(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        syn_measure.(glacier)   = diag(sampledtemp);
        
        real_measure_tmp.(glacier) = real_measure.(glacier)(~isnan(syn_measure.(glacier)));

            syn_measure_tmp.(glacier) = syn_measure.(glacier)(~isnan(syn_measure.(glacier)));
        rmseBMA.(namesP{p}).(glacier)(ss,mc) = sqrt(mean((syn_measure_tmp.(glacier)-real_measure_tmp.(glacier)).^2));
        end
        
        
        % BASIC LR        
        for g = 3;        glacier = char(options.glacier(g));

        swe	    = WBinputN.(glacier)(:,1);
        Xt      = struct2array(TOPOinput.(namesP{p}).(glacier));
        X       = [ones(length(Xt),1), Xt];

        coeffs = regress(swe, X);
        coeffsLR(ss).(namesP{p})(mc).(glacier)   = coeffs;

        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;
        swe_pred = swe_pred.*AblationArea.(glacier); 
        predLR(ss).(namesP{p})(mc).(glacier) = swe_pred;
        
        sampledtemp     = swe_pred(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        syn_measure     = diag(sampledtemp);
        
        real_measure_tmp = real_measure.(glacier)(~isnan(syn_measure));
        syn_measure_tmp = syn_measure(~isnan(syn_measure));
        
        rmseLR.(namesP{p}).(glacier)(ss,mc) = sqrt(mean((syn_measure_tmp-real_measure_tmp).^2));
        end
end
end
end

%% Plot coeffs

pattern = 'CentreTransect';
m = 1;
for i=1:10
    for g = 3;        glacier = char(options.glacier(g));
     subplot(5,2,m)
%     plot_coeffs = [coeffSYN(i,:)',coeffsLR(nn).(pattern)(i).(glacier)(2:end), coeffsBMA(nn).(pattern)(i).(glacier)(1:7);
%                     intercept(g), coeffsLR(nn).(pattern)(i).(glacier)(1), coeffsBMA(nn).(pattern)(i).(glacier)(8)];
    plot_coeffs = [coeffSYN(i,:)',coeffsLR(nn).(pattern)(i).(glacier)(2:end);
                    intercept(g), coeffsLR(nn).(pattern)(i).(glacier)(1)];
   bar(plot_coeffs)
    params={'z','\alpha','m','N','\kappa','Sx','intercept'};
    set(gca,'xticklabel',params)
    ylabel('\beta')
    title(['Model ',num2str(i),' ',glacier])
    m=m+1;
%     legend('Model','Simple LR','BMA')
    end
end

%% Plot estimated distributions

pattern = 'CentreTransect';

for m =10%1:num_models
    figure(1);
    PlotTopoParameter(predLR(nn).(pattern)(m), 'none', 'WB (m w.e.)', pUTM.(pattern), 'black', 'none')
    figure(2);
    PlotTopoParameter(predBMA(nn).(pattern)(m), 'none', 'WB (m w.e.)', pUTM.(pattern), 'black', 'none')
    
end

%% 

pattern = 'CentreTransect';
m = 1;
for i=1:5
    for g = 3;        glacier = char(options.glacier(g));
    subplot(5,2,m)
    actual = snowdist_model(i).(glacier).*AblationArea.(glacier);
    imagesc(actual); colorbar
        caxis([0 2]); caxis(caxis);
    title({['Model ',num2str(i)],['Model ', glacier, ' B_W=', num2str(nanmean(actual(:)))]})
    axis('off')
    
    subplot(5,2,m+1)
    imagesc(predLR(nn).(pattern)(i).(glacier)); colorbar
        caxis([0 2]); caxis(caxis);
    title({['Model ',num2str(i)],['Simple LR ', glacier, ' B_W=', num2str(nanmean(predLR(nn).(pattern)(i).(glacier)(:)))]})
    axis('off')
    
%     subplot(5,3,m+2)
%     imagesc(predBMA(nn).(pattern)(i).(glacier)); colorbar
%         caxis([0 2]); caxis(caxis);
%     title({['Model ',num2str(i)],['BMA ', glacier, ' B_W=', num2str(nanmean(predBMA(nn).(pattern)(i).(glacier)(:)))]})
%     axis('off')
    
    m=m+2;
%     legend('Simple LR','BMA')
    end
end


%%
glacier='G13';
for i=1:10
x(i)=nanmean(nanmean(snowdist_model(i).(glacier).*AblationArea.(glacier)));
end


%% Plot beta distributions

elev_range  = [0.585, 0.810];
dc_range    = [0, 0.015];
aspect_range = [0, 0.093];
slope_range = [-0.09, 0.277];
N_range     = [-0.414, 0.172];
curv_range  = [-0.077, 0.02];
sx_range    = [-0.294, 0.260];

intercept   = [0.6205, 0.2627, 0.2354];

% Normal random beta values

x = [-1:.01:1.5];
elev_beta   = normpdf(x, mean(elev_range), std(elev_range));
dc_beta     = normpdf(x, mean(dc_range), std(dc_range));
aspect_beta = normpdf(x, mean(aspect_range), std(aspect_range));
slope_beta  = normpdf(x, mean(slope_range), std(slope_range));
N_beta      = normpdf(x, mean(N_range),std(N_range));
curv_beta   = normpdf(x, mean(curv_range), std(curv_range));
sx_beta     = normpdf(x, mean(sx_range),   std(sx_range));

plot(x,elev_beta); hold on
plot(x,dc_beta); hold on
plot(x,aspect_beta); hold on
plot(x,slope_beta); hold on
plot(x,N_beta); hold on
plot(x,curv_beta); hold on
plot(x,sx_beta); hold on
legend('z','d_c','\alpha','m','N','\kappa','Sx');
title('\beta distirbutions')
xlabel('\beta')
ylabel('Probability')


%% Number of synthetic models
load PaperII_LR_runtest

pattern = 'CentreTransect';
glacier = 'G13';
ss = 18;
edges = 0:0.03:0.6;

histogram(rmseLR500.(pattern).(glacier)(ss,:), edges); hold on
histogram(rmseLR200.(pattern).(glacier)(ss,:), edges); hold on
histogram(rmseLR100.(pattern).(glacier)(ss,:), edges); hold on
legend('100','200','500')
title('Number of synthetic models')
ylabel('Frequency')
xlabel('RMSE')

display([mean(rmseLR500.(pattern).(glacier)(ss,:)), std(rmseLR500.(pattern).(glacier)(ss,:))])
display([mean(rmseLR200.(pattern).(glacier)(ss,:)), std(rmseLR200.(pattern).(glacier)(ss,:))])
display([mean(rmseLR100.(pattern).(glacier)(ss,:)), std(rmseLR100.(pattern).(glacier)(ss,:))])