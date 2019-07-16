clear
load PaperII_Syntheic.mat
load PaperII_RegularSampling
load PaperII_AblationArea.mat AblationArea AccumulationArea

load PaperII_SynSnowDistModel_Pul_Abl_DDA

% file_path = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/';
file_path = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/';

sampling_number = 200;
num_models = 1:200;
nn = 6;

only_ablation = true;

file_name_save = 'PII_Syn_July5_Pul_abl.mat';

Bw_abl = [0.59, 0.34, 0.27]; % ablation area only
Bw_full = [0.59, 0.57, 0.38]; % whole glacier
if only_ablation
    Bw_obs = Bw_abl;
else
    Bw_obs = Bw_full;
end
    p_acc = [0.36,    0.3,        0.308];
    p_abl = [0.3422,  0.353,     0.3692];

%% SAMPLING THEORETICAL FIELD FROM PATTERNS

 %Make matrix with utm of each grid cell
for g = 1:3;    glacier = options.glacier{g};
minE = min(options.rig.(glacier)(:,1));
minN = min(options.rig.(glacier)(:,2));

    nE = options.mapsize(g,2);            
    nN = options.mapsize(g,1);
    
    utmGridE.(glacier) = repmat([1:nE]*40+minE,nN,1);   
    utmGridN.(glacier) = repmat([nN:-1:1]'*40+minN,1,nE); 
end

for mc = num_models%1:num_models

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
namesP = fieldnames(pWB);
% namesP = {'Centreline'};

real_measure = SampledCell(snowdist_model(mc));

for p = 1:length(namesP)    

for ss = nn
    
    for g = 1:3;        glacier = char(options.glacier(g));
       display([' Sample size: ',num2str(ss),' Pattern: ',namesP{p}, ' Run:',num2str(mc), ' Glacier ', num2str(g)])

% DETERMINISTIC
        sampleN = 1;
    [WBinput, TOPOinput_temp] = deterministic_locations([pWB.(namesP{p}).(glacier), pUTM.(namesP{p}).(glacier)],...
                                             pTOPO.(namesP{p}).(glacier), ss);
        subUTM.(namesP{p}).(glacier) = WBinput(:,2:3);
     
        ff = fieldnames(TOPOinput_temp);
        TOPOinput = zeros(ss,length(ff));
    for i = 1:length(ff);    fname = ff{i};
    TOPOinput(:,i)  = TOPOinput_temp.(fname); end
    clear syn_measure_tmp syn_measure real_measure_tmp
        swe	    = WBinput(:,1);
            swe     = swe + normrnd( 0, options.zzstd(g), size(swe,1),1);
            swe(swe<0) = 0;
        Xt      = TOPOinput;
        X       = [ones(size(Xt,1),1), Xt];

        coeffs = regress(swe, X);
        coeffsLR_det(ss).(namesP{p})(mc).(glacier)   = coeffs;

        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;
        
        swe_pred(AccumulationArea.(glacier)==1) = swe_pred(AccumulationArea.(glacier)==1)*p_acc(g)./p_abl(g);
%         synBw_det_wAccum.(namesP{p}).(glacier)(ss,mc,sampleN) = nanmean(swe_pred(:));

        if only_ablation
            swe_pred = swe_pred.*AblationArea.(glacier); 
        end
        synBw_det.(namesP{p}).(glacier)(ss,mc,sampleN) = nanmean(swe_pred(:));
        
        % RMSE
        real_snow = snowdist_model(mc).(glacier)(:);
        syn_snow = swe_pred(:);
        
        xNaN = ~any(isnan([real_snow,syn_snow]),2);
        
        R = corrcoef(real_snow(xNaN), syn_snow(xNaN));
        R2_det.(namesP{p}).(glacier)(ss,mc,sampleN) = R(2,1)^2;
        
        synRMSE_det.(namesP{p}).(glacier)(ss,mc,sampleN) = ...
            sqrt(nanmean(nanmean((real_snow(xNaN)-syn_snow(xNaN)).^2)));



%RANDOM UNIFORM
    for sampleN = 1:sampling_number

[WBinput,TOPOinput_temp] = random_uniform([pWB.(namesP{p}).(glacier),fullUTM.(namesP{p}).(glacier)],...
                                             pTOPO.(namesP{p}).(glacier), ss); 
     subUTM.(namesP{p}).(glacier) = WBinput(:,2:3);
     
        ff = fieldnames(TOPOinput_temp);
        TOPOinput = zeros(ss,length(ff));
    for i = 1:length(ff);    fname = ff{i};
    TOPOinput(:,i)  = TOPOinput_temp.(fname); end

    clear syn_measure_tmp syn_measure real_measure_tmp
    
        % BASIC LR        

        swe	    = WBinput(:,1);
            swe     = swe + normrnd( 0, options.zzstd(g), size(swe,1),1);
            swe(swe<0) = 0;
        Xt      = TOPOinput;
        X       = [ones(size(Xt,1),1), Xt];
        
        coeffs = regress(swe, X);
        coeffsLR_rand(ss).(namesP{p})(mc).(glacier)   = coeffs;
        
%         display([swe,X])
%         display(coeffs)
%         pause

        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;
        
        swe_pred(AccumulationArea.(glacier)==1) = swe_pred(AccumulationArea.(glacier)==1)*p_acc(g)./p_abl(g);
%         synBw_rand_wAccum.(namesP{p}).(glacier)(ss,mc,sampleN) = nanmean(swe_pred(:));

        if only_ablation
            swe_pred = swe_pred.*AblationArea.(glacier); 
        end
        predLR(ss).(namesP{p})(mc).(glacier) = swe_pred;
        synBw_rand.(namesP{p}).(glacier)(ss,mc,sampleN) = nanmean(swe_pred(:));
        
        % RMSE
        real_snow = snowdist_model(mc).(glacier)(:);
        syn_snow = swe_pred(:);
        
        xNaN = ~any(isnan([real_snow,syn_snow]),2);
        
        R = corrcoef(real_snow(xNaN), syn_snow(xNaN));
        R2_rand.(namesP{p}).(glacier)(ss,mc,sampleN) = R(2,1)^2;
        
        synRMSE_rand.(namesP{p}).(glacier)(ss,mc,sampleN) = ...
            sqrt(nanmean(nanmean((real_snow(xNaN)-syn_snow(xNaN)).^2)));
    end

    end
        
end
end
end


%% Calculate best RMSE using all data

% load PaperII_FinalLRruns
for mc = num_models

%Additing noise to distributed WB
real_measure = SampledCell(snowdist_model(mc));

topoparam = fieldnames(topo_full.G4);
for num_param = 1:length(topoparam)
    for g = 1:3;        glacier = char(options.glacier(g));
    tmp_topo.(glacier) = topo_full.(glacier).(topoparam{num_param});
    end    
    tmp_sampled.(topoparam{num_param}) = SampledCell(tmp_topo);
end
for g = 1:3;        glacier = char(options.glacier(g));
    for num_param = 1:length(topoparam)
sampled_topo.(glacier).(topoparam{num_param}) = tmp_sampled.(topoparam{num_param}).(glacier);
    end
end        
        % BASIC LR        
for g = 1:3;        glacier = char(options.glacier(g));

        swe	    = real_measure.(glacier)(:,1);
            swe     = swe + normrnd( 0, options.zzstd(g), size(swe,1),1);
            swe(swe<0) = 0;
        Xt      = struct2array(sampled_topo.(glacier));
        X       = [ones(length(Xt),1), Xt];

        coeffs = regress(swe, X);
        tmp_best_coeffsLR.(glacier)(mc,:)   = coeffs;

        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for num_param = 1:length(betaCoeff)
            param      = topoCoeff{num_param};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_param);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;
        if only_ablation
            swe_pred = swe_pred.*AblationArea.(glacier);
        end
        best_predLR.(glacier) = swe_pred;
        
        % RMSE
        real_snow = snowdist_model(mc).(glacier)(:);
        syn_snow = swe_pred(:);
        
        xNaN = ~any(isnan([real_snow,syn_snow]),2);
        
        R = corrcoef(real_snow(xNaN), syn_snow(xNaN));
        tmp_best_R2.(glacier)(mc) = R(2,1)^2;
        
%         tmp_best_rmseLR.(glacier)(mc) = sqrt(mean((real_snow(xNaN)-syn_snow(xNaN)).^2));
        tmp_best_rmseLR.(glacier)(mc) = sqrt(nanmean(nanmean((snowdist_model(mc).(glacier)-swe_pred).^2)));
        
%         figure(1); clf
%         subplot(1,3,1)
%         h=imagesc(snowdist_model(mc).(glacier));
%         set(h,'alphadata',~isnan(snowdist_model(mc).(glacier))); caxis([0 1]); caxis(caxis); axis off 
%         
%         subplot(1,3,2)
%         h=imagesc(swe_pred);
%         set(h,'alphadata',~isnan(swe_pred)); caxis([0 1]); caxis(caxis); axis off 
%         
%         subplot(1,3,3)
%         h=imagesc(abs(snowdist_model(mc).(glacier)-swe_pred));
%         set(h,'alphadata',~isnan(snowdist_model(mc).(glacier)-swe_pred)); caxis([0 0.5]); caxis(caxis); axis off 
%         
%         display(tmp_best_rmseLR.(glacier)(mc))
%         pause
end   
end

for g = 1:3;        glacier = char(options.glacier(g));
best_rmseLR.(glacier) = mean(tmp_best_rmseLR.(glacier));
best_coeffsLR.(glacier) = mean(tmp_best_coeffsLR.(glacier));
best_R2.(glacier) = mean(tmp_best_R2.(glacier));
end

%%
save(file_name_save,'synBw_rand','synBw_det','R2_rand','R2_det',...
    'synRMSE_rand','synRMSE_det','best_rmseLR','best_coeffsLR','best_R2')