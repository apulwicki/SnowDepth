%% THEORETICAL - All n for WB

file_path = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/';
% file_path = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/';

%% Get synthetic snow distributions
load PaperII_Syntheic.mat

% Beta coeff ranges
elev_range  = [0.05, 0.810];
sx_range    = [-0.294, 0.260];
curv_range  = [-0.077, 0.02];
slope_range = [-0.09, 0.277];
intercept   = [0.6205, 0.2627, 0.2354];

% Normal random beta values
num_syn_models = 100;

elev_beta   = normrnd(mean(elev_range), std(elev_range), [1,num_syn_models]);
sx_beta     = normrnd(mean(sx_range),   std(sx_range),   [1,num_syn_models]);
curv_beta   = normrnd(mean(curv_range), std(curv_range), [1,num_syn_models]);
slope_beta  = normrnd(mean(slope_range),std(slope_range),[1,num_syn_models]);

    % Remove dc, aspect and Northness
    for g = 1:3;    glacier = options.glacier{g};
    topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
    end
    
for g=1:3; glacier = options.glacier{g};
AblationArea.(glacier)(AblationArea.(glacier)==-0.1)=NaN;
AblationArea.(glacier)(~isnan(AblationArea.(glacier)))=1;
end


for m = 1:num_syn_models

for g = 1:3;    glacier = options.glacier{g};
    syn_model   = repmat(intercept(g), options.mapsize(g,:));
    betaCoeff   = [elev_beta(m), sx_beta(m), curv_beta(m), slope_beta(m)];   
    topoCoeff   = fieldnames(topo_full.G4);
    
    for c = 1:length(betaCoeff)
        param     = topoCoeff{c};
        sweT      = topo_full.(glacier).(param)*betaCoeff(c);
        syn_model = syn_model + sweT;
    end
        
    syn_model(syn_model<0) = 0; %Set min to 0
    syn_model = syn_model.*AblationArea.(glacier);
    
    snowdist_model(m).(glacier) = syn_model;
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

for mc = 1:num_syn_models

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

    
%%
namesP = fieldnames(pWB);
%  namesP = {'Circle'};

real_measure = SampledCell(snowdist_model(mc));

for p = 1:length(namesP)

for ss = 6:55
   display([' Sample size: ',num2str(ss),' Pattern: ',namesP{p}, ' Run:',num2str(mc)])

   [WBinput, TOPOinput, UTMinput] = SubsetSampleSize( pWB, pTOPO, pUTM, ss );


    %Add some noise
    WBinputN = WBnoise(WBinput.(namesP{p}),'low');
    
        for g = 1:3;        glacier = char(options.glacier(g));
        
        % BASIC LR
        swe	    = WBinputN.(glacier)(:,1);
        Xt      = struct2array(TOPOinput.(namesP{p}).(glacier));
        X       = [ones(length(Xt),1), Xt];

        % Get coefficients
        coeffs = regress(swe, X);
%         synMLR(ss).(namesP{p})(mc).(glacier)   = coeffs;
        
        %Predict
        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
            %multiply coeffs and add them
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
            %Set min to 0
        swe_pred(swe_pred<0) = 0;
        swe_pred = swe_pred.*AblationArea.(glacier); 
        
%         SynPred(ss).(namesP{p})(mc).(glacier) = swe_pred;

        % Calculate RMSE
        true_dist = snowdist_model(mc).(glacier);
        pred_dist = swe_pred;
        
        
        synRMSE.(namesP{p}).(glacier)(ss,mc) = sqrt(nanmean(nanmean((pred_dist-true_dist).^2)));
        end

end
end
end

% synRMSE_high = synRMSE;
% synCOEF_high = synMLR;
test = test2;
synPred_utm = UTMinput;
save('PII_FastRuns.mat','SynPred', 'synPred_utm', 'snowdist_model', '-append')
% save('PII_FastRuns.mat','synRMSE_high','synCOEF_high','-append')

%% Best RMSE for all synthetic models

for g = 1:3;        glacier = char(options.glacier(g));
    elev.(glacier)   = topo_full.(glacier).elevation;
    slope.(glacier)  = topo_full.(glacier).slope;
    curve.(glacier)  = topo_full.(glacier).curvature;
    Sx.(glacier)     = topo_full.(glacier).Sx;
end

% [ ~, measure_topo, ~ ]    = ObsInCell(fullSWE.S2.input, topo_sampled);

bestRMSE = zeros(num_syn_models,3);
% [WBinput, TOPOinput, UTMinput] = SubsetSampleSize( pWB, pTOPO, pUTM, 200);

for mc = 1:num_syn_models
    %Measurement locations
    measure_points  = SampledCell(snowdist_model(mc));
%     measure_points  = WBnoise(measure_points,'low');    
%     measure_HC = WBnoise(WBinput.HourCircle,'low');
%     measure_CT = WBnoise(WBinput.CentreTransect,'low');
    
    measure_elev    = SampledCell(elev);
    measure_slope   = SampledCell(slope);
    measure_curve   = SampledCell(curve);
    measure_Sx      = SampledCell(Sx);
    
    for g = 1:3;        glacier = char(options.glacier(g));
        
%         measure_points = [measure_HC.(glacier); measure_CT.(glacier)];

%         elev    = [TOPOinput.HourCircle.(glacier).elevation; TOPOinput.CentreTransect.(glacier).elevation];
%         slope   = [TOPOinput.HourCircle.(glacier).slope; TOPOinput.CentreTransect.(glacier).slope];
%         curve   = [TOPOinput.HourCircle.(glacier).curvature; TOPOinput.CentreTransect.(glacier).curvature];
%         Sx      = [TOPOinput.HourCircle.(glacier).Sx; TOPOinput.CentreTransect.(glacier).Sx];
        
        % BASIC LR
        swe	    = measure_points.(glacier);
%         Xt      = [measure_topo.(glacier).elevation measure_topo.(glacier).slope...
%                    measure_topo.(glacier).curvature measure_topo.(glacier).Sx];
        Xt      = [measure_elev.(glacier) measure_slope.(glacier) measure_curve.(glacier) measure_Sx.(glacier)];
%         Xt      = [elev, slope, curve, Sx];
        X       = [ones(length(Xt),1), Xt];

        coeffs = regress(swe, X);        
        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for n = 1:length(betaCoeff)
            sweT       = topo_full.(glacier).(topoCoeff{n})*betaCoeff(n);
            swe_pred   = swe_pred + sweT;
        end
            %Set min to 0
        swe_pred(swe_pred<0) = 0;
        swe_pred = swe_pred.*AblationArea.(glacier); 
        
        % Calculate RMSE
        true_dist = snowdist_model(mc).(glacier);
        pred_dist = swe_pred;
        
        disp(sqrt(nanmean(nanmean((pred_dist-true_dist).^2))))
        bestRMSE(mc,g) = sqrt(nanmean(nanmean((pred_dist-true_dist).^2)));
    end
end
%% PLOT

load Full.mat  fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS
load Patterns.mat pUTM
load PaperII_AblationArea.mat
load PII_FastRuns.mat

    namesP = fieldnames(synRMSE);  order = [2 3 1 4 5 6]; namesP = namesP(order);
    N1 = {'Centreline',''}; N2 = {'Centreline &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Random',''}; 
    namesPfull = [N1; N2; N3; N4;N5;N6];
    C =[     0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.0588    0.3490    0.1216];
    realGrid    = ObsInCell(fullSWE.S2.input, topo_sampled);
     ConvTable   = nan(3,6);  ConvTable = array2table(ConvTable,'VariableNames',namesP);  
     VarTable    = nan(3,6);  VarTable  = array2table(VarTable,'VariableNames',namesP);


    numPoints = 6:45;
    RMSEfull = mean(bestRMSE);
    ela_ind = [1 7; 8 15; 16 23];
    ela_m = [0 0; -0.04*10^5 -0.0088*10^6; 0.04*10^5 0.0086*10^6];
    
  % Figure  
     clf; n = 1;
[ha, ~] = tight_subplot(3,length(namesP),[0.01 .01],[.08 0],[.05 0]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
    
meanWB      = mean(synRMSE_low.(namesP{p}).(glacier)(numPoints,:),2);
stdWB_low   = std(synRMSE_low.(namesP{p}).(glacier)(numPoints,:),[],2);
stdWB_high  = std(synRMSE_high.(namesP{p}).(glacier)(numPoints,:),[],2);

axes(ha(n))
plot(numPoints,meanWB,'LineWidth',0.5,'Color',C(p,:)); hold on
    upper_low = meanWB + stdWB_low;
    lower_low = meanWB - stdWB_low;
fill([numPoints flip(numPoints)],[upper_low',flip(lower_low')],...
     C(p,:),'FaceAlpha',0.4,'EdgeColor','none'); hold on
    upper_high = meanWB + stdWB_high;
    lower_high = meanWB - stdWB_high;
fill([numPoints flip(numPoints)],[upper_high',flip(lower_high')],...
     C(p,:),'FaceAlpha',0.2,'EdgeColor','none'); hold on

 
    smoothSize = 3;
    %Smooth mean  
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = meanWB(h:end-(smoothSize-h));  end
    meanSmooth = mean(M);    
    plot(numPoints((h-1)/2:end-(h-1)/2-1), meanSmooth,'LineWidth',1.8,'Color',C(p,:)); hold on
    
    %Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%
%      estGrid = SampledCell(sweBMS_alldata);
%      RMSEfull = sqrt(mean((estGrid.(glacier)-realGrid.(glacier)(:,1)).^2));
    plot([min(numPoints) max(numPoints)],[RMSEfull(g),RMSEfull(g)],'--k')

    % 5% of GW mean
%     good = find(T-RMSEfull<0.05,1);  
    good = find(meanSmooth-RMSEfull(g)<0.05,1);
%       good = find(2*stdWB_low<0.05,1);
        if ~isempty(good)
        good = good + min(numPoints);%-1+(smoothSize-1)/2;
        plot([good good],[0 2],'k-.','LineWidth',0.05)
        ConvTable{g,p}  = good;
        end
%     %High STD within 25% of mean
%         M = zeros(smoothSize, length(meanWB)-smoothSize+1);
%     for h = 1:smoothSize;  M(h,:) = stdWB(h:end-(smoothSize-h));  end
%     TN = mean(M);    
%     N10 = find(TN/RMSEfull<0.25,1); 
%         if ~isempty(N10) 
%         N10 = N10+min(numPoints)-1+(smoothSize-1)/2;
%         plot([N10 N10],[0 1.2],'-.k','LineWidth',0.75)
%         VarTable{g,p}   = N10;
%         end

    %titles
    if g==1; title(namesPfull(p,:)); end
    %y labels
    if g ==1 && p==1; ylabel('G4 RMSE (m w.e.)'); set(gca,'XTickLabel',[]);
    elseif g ==2 && p==1; ylabel('G2 RMSE (m w.e.)'); set(gca,'XTickLabel',[]);
    elseif g ==3 && p==1; ylabel('G13 RMSE (m w.e.)'); 
    else set(gca,'YTickLabel',[]);
    end
    if n==15; xlabel('                            Sample size'); end
    if g==1 || g == 2; set(gca,'XTickLabel',[]); end
ylim([0 1])
xlim([6 45])
n = n+1;

ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);

    if n==2; yInd = 0.002; else; yInd = 0.14; end
%     yInd = 0.14;
    if g ==1; sizeG = 0.12; xoff = 0.04; yInd = yInd+0.02; else sizeG = 0.14; xoff = 0.036; end
axes('position',[ax(1)+xoff ax(2)+yInd sizeG sizeG])
pInd = 1:10:length(pUTM.(namesP{p}).(glacier));
if g == 3; fill(options.rig.(glacier)(1:304,1),options.rig.(glacier)(1:304,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
else; fill(options.rig.(glacier)(:,1),options.rig.(glacier)(:,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
end
plot(pUTM.(namesP{p}).(glacier)(pInd,1),pUTM.(namesP{p}).(glacier)(pInd,2),'k.','MarkerSize',1.5); hold on
plot(options.rig.ELA(ela_ind(g,1):ela_ind(g,2),1)+ela_m(g,1), options.rig.ELA(ela_ind(g,1):ela_ind(g,2),2)+ela_m(g,2),'-.k');
    axis off; axis equal
end
end 

%     saveFIG_HP('PII_AA_SyntheticData',2,12)


%% PLOT - DIFFERENCE MAP RMSE

% load PII_FastRuns.mat SynPred snowdist_model synPred_utm
load Full.mat fullLR
run OPTIONS
load PaperII_AblationArea.mat AblationArea

    namesP = fieldnames(SynPred);  order = [2 3 1 4 5 6]; namesP = namesP(order);
    N1 = {'Centreline',''}; N2 = {'Centreline &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Random',''};
    namesPfull = [N1; N2; N3; N4;N5;N6];


figure(1); clf
SizeToView = 20;
n = 1;
nRuns = 100;%num_syn_models;

[ha, pos] = tight_subplot(3,length(namesP),[0 0],[0 0],[0 0]);
annStart = [0.061 0.77 .1 .1];

for p = 1:length(namesP)
for g = 1:3; glacier = options.glacier{g};
for mc = 1:nRuns
    pred_Syn_full.(namesP{p}).(glacier)(:,:,mc) = SynPred(SizeToView).(namesP{p})(mc).(glacier);     
end  
end
end


for g = 1:3 
for p = 1:length(namesP)
    glacier = options.glacier{g};
%Calculate relative uncertainty and plot
    C = snowdist_model(mc).(glacier);    C = repmat(C,[1 1 nRuns]);
    D = pred_Syn_full.(namesP{p}).(glacier)-C;
    RMSEdist = sqrt(sum(D.^2,3)/nRuns);     
    
    RMSEdist(AblationArea.(glacier)==-0.1) = -0.1;
    
        c = nanmean(snowdist_model(mc).(glacier)(:)); 
        d = nanmean(pred_Syn_full.(namesP{p}).(glacier),2);
        d = nanmean(d,1);   d = d(:)-c;
    RMSEglacier = sqrt(sum(d.^2)/nRuns);
    
axes(ha(n));    
h(n) = imagesc(RMSEdist); hold on
    c = cbrewer('seq', 'BuPu', 100,'PCHIP');
    c(1,:) = [199, 201, 204]/255;
    colormap(c);         set(h(n),'alphadata',~isnan(RMSEdist));
    axis square; axis off;
    caxis([-0.1 1.5])
    
%Plot sampling locations
    E = (synPred_utm.(namesP{p}).(glacier)(:,1)-min(options.rig.(glacier)(:,1)))/40;
        minN = min(options.rig.(glacier)(:,2));
        Ng = (options.rig.(glacier)(:,2) - minN)/40; 
    N = (synPred_utm.(namesP{p}).(glacier)(:,2)-minN)/40; N = max(Ng)-N;
plot(E,N,'k.','MarkerSize',4.5)

%Add RMSE of glacier-wide WB
    annotation('textbox',[annStart(1)+(p-1)*0.17 annStart(2)-(g-1)*0.3 .1 .1],...
        'String', [num2str(round(RMSEglacier,2),'%1.2f'),' m w.e.'],'EdgeColor','none','FontWeight','bold')

    if g==1; title(namesPfull(p,:)); end

    n = n+1;
end
end

% 	saveFIG_HP('PII_AA_DistributedRMSE',2,12)


%% nc calculation

load PII_FastRuns.mat
load PaperII_AblationArea.mat sweBMS_alldata
load Full.mat fullLR fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS

numPoints   = 6:57;
namesP      = fieldnames(synRMSE_low);
nc_table    = nan(3,6);     nc_table = array2table(nc_table,'VariableNames',namesP);  
nv_table    = nan(3,6);     nv_table = array2table(nv_table,'VariableNames',namesP);  
realGrid    = ObsInCell(fullSWE.S2.input, topo_sampled);
RMSEfull    = mean(bestSynRMSE);


for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
    
meanWB  = mean(synRMSE_low.(namesP{p}).(glacier)(numPoints,:),2);
stdWB   = std(synRMSE_low.(namesP{p}).(glacier)(numPoints,:),[],2);

    %Smooth mean 
    smoothSize = 3;
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = meanWB(h:end-(smoothSize-h));  end
    meanSmooth = mean(M); 
    %Smooth std  
        M = zeros(smoothSize, length(stdWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = stdWB(h:end-(smoothSize-h));  end
    stdSmooth = mean(M);    



nc = find(meanSmooth-RMSEfull(g)<RMSEfull(g)*0.1,1);
if isempty(nc); nc = nan; 
else; nc = nc + min(numPoints); end
nc_table{g,p}   = nc;

nv = find(stdSmooth+meanSmooth-RMSEfull(g)<RMSEfull(g)*0.75,1);
if isempty(nv); nv=nan; 
else; nv = nv + min(numPoints); end
nv_table{g,p}   = nv;
    
end 
end