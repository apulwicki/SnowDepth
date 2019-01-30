%% Figure 2 - Synthetic sampling designs 
    load PaperII_AblationArea.mat AblationArea
    load Patterns.mat pUTM 
    load Full.mat fullLR
    run OPTIONS
NumSubPoints = 15;

% Blank our accum area
for g = 1:3;    glacier = options.glacier{g};
    abla_area.(glacier) = fullLR.S2.(glacier);
    abla_area.(glacier)(AblationArea.(glacier)==-0.1)=-0.1;
end


    P = fieldnames(pUTM);
for t = 1:length(P)
    clear pattern*
for g = 1:3;    glacier = options.glacier{g};
        if t == 6; pUTM.(P{t}).(glacier)(:,1) = []; end
patternFULL.(glacier)(:,2:3) = pUTM.(P{t}).(glacier);

I = floor(linspace(1,size(pUTM.(P{t}).(glacier),1),NumSubPoints));
patternSUB.(glacier)(:,2:3) = pUTM.(P{t}).(glacier)(I,:);
end

figure(1)
PlotTopoParameter_DD(abla_area, 'WB (m w.e.)', patternFULL, patternSUB, 'black')
	saveFIG_HP(['SampleDesign_',P{t}],1,8.6)
end

%% Figure 3 - Synthetic data WB 
%    clear
load Patterns.mat 
load Patterns.mat pUTM
load Full.mat fullLR 
% load PaperII_FinalLRruns rmseLR mpeLR
run OPTIONS

    namesP = fieldnames(SynObs_High);  order = [2 3 1 4 5 7]; namesP = namesP(order);
    N1 = {'Centreline',''}; N2 = {'Centreline &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Random',''};
    namesPfull = [N1; N2; N3; N4;N5;N6];
    pUTM.Random.G2(:,1) = []; pUTM.Random.G4(:,1) = []; pUTM.Random.G13(:,1) = [];   
    C =[     0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.0588    0.3490    0.1216];
    ela_ind = [1 7; 8 15; 16 23];
    ela_m = [0 0; -0.04*10^5 -0.0088*10^6; 0.04*10^5 0.0086*10^6];
    
  % Figure  
     clf; n = 1;
     ConvTable   = nan(3,6);  ConvTable = array2table(ConvTable,'VariableNames',namesP);  
     VarTable    = nan(3,6);  VarTable  = array2table(VarTable,'VariableNames',namesP);

[ha, ~] = tight_subplot(3,length(namesP),[0.01 .01],[.08 0],[.05 0]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
%    numPoints = 8:size(SynObs_High.(namesP{p}).(glacier),1);
   numPoints = 5:30;

%     M1   = mean(rmseLR_nonoise.(namesP{p}).(glacier)(numPoints,:),2);
%     M2   = mean(rmseLR_low.(namesP{p}).(glacier)(numPoints,:),2);   
% meanWB   = mean([M1 M2],2);
meanWB = mean(rmseLR.(namesP{p}).(glacier)(numPoints,:),2);
stdWB  = std(rmseLR.(namesP{p}).(glacier)(numPoints,:),[],2);
% stdWB_L  = std(rmseLR_low.(namesP{p}).(glacier)(numPoints,:),[],2);


axes(ha(n));    
plot(numPoints,meanWB,'LineWidth',0.5,'Color',C(p,:)); hold on
    upper = meanWB + stdWB;
    lower = meanWB - stdWB;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.5,'EdgeColor','none'); hold on
%      upper = meanWB + stdWB_L;
%     lower = meanWB - stdWB_L;
% fill([numPoints flip(numPoints)],[upper',flip(lower')],...
%      C(p,:),'FaceAlpha',0.4,'EdgeColor','none'); hold on

%Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%
   
    GWmean = best_rmseLR.(glacier); %nanmean(fullLR.S2.(glacier)(:));
    
    smoothSize = 7;
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = meanWB(h:end-(smoothSize-h));  end
    T = mean(M);    
    
    plot(numPoints((h-1)/2:end-(h-1)/2-1), T,'LineWidth',1.6,'Color',C(p,:)); 
    plot([min(numPoints) max(numPoints)],[GWmean,GWmean],'--k')

    % 5% of GW mean
    good = find(T-GWmean<0.05,1);   
        if ~isempty(good)
        good = good + min(numPoints)-1+(smoothSize-1)/2;
%         plot([good good],[0 2],'k:','LineWidth',0.05)
        ConvTable{g,p}  = good;
        end
    %High STD within 25% of mean
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = stdWB(h:end-(smoothSize-h));  end
    TN = mean(M);    
    N10 = find(TN/GWmean<0.75,1); 
        if ~isempty(N10) 
        N10 = N10+min(numPoints)-1+(smoothSize-1)/2;
%         plot([N10 N10],[0 1.2],'-.k','LineWidth',0.75)
        VarTable{g,p}   = N10;
        end
        

    %titles
    if g==1; title(namesPfull(p,:)); end
    %y labels
    if g ==1 && p==1; ylabel('G4 MPE (%)'); set(gca,'XTickLabel',[]);
    elseif g ==2 && p==1; ylabel('G2 MPE (%)'); set(gca,'XTickLabel',[]);
    elseif g ==3 && p==1; ylabel('G13 MPE (%)'); 
    else set(gca,'YTickLabel',[]);
    end
    if n==15; xlabel('                            Sample size'); end
    if g==1 || g == 2; set(gca,'XTickLabel',[]); end
ylim([0 0.3])
xlim([0 30])
n = n+1;

ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);

    if n==223; yInd = 0.002; else yInd = 0.14; end
    if g ==1; sizeG = 0.12; xoff = 0.04; yInd = yInd+0.02; else sizeG = 0.14; xoff = 0.036; end
axes('position',[ax(1)+xoff ax(2)+yInd sizeG sizeG])
pInd = 1:length(pUTM.(namesP{p}).(glacier));
if g == 3; fill(options.rig.(glacier)(1:304,1),options.rig.(glacier)(1:304,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
else fill(options.rig.(glacier)(:,1),options.rig.(glacier)(:,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
end
plot(pUTM.(namesP{p}).(glacier)(pInd,1),pUTM.(namesP{p}).(glacier)(pInd,2),'k.','MarkerSize',1.5);
plot(options.rig.ELA(ela_ind(g,1):ela_ind(g,2),1)+ela_m(g,1), options.rig.ELA(ela_ind(g,1):ela_ind(g,2),2)+ela_m(g,2),'-k');
    axis off; axis equal
end
end 

	%saveFIG_HP('SyntheticObsWB',2,12)
%     saveFIG_HP('Pulwicki_Fig3_MPE',2,12)

%% Figure 5 - Real data
    clear
load Patterns.mat 
load Full.mat fullLR fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS

    namesP = fieldnames(DataObs_RMSE);  order = [2 3 1 4 5 6]; namesP = namesP(order);
    N1 = {'Midline',''}; N2 = {'Mid &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Random',''}; 
    namesPfull = [N1; N2; N3; N4;N5;N6];
    pUTM.Random.G2(:,1) = []; pUTM.Random.G4(:,1) = []; pUTM.Random.G13(:,1) = [];   
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



  % Figure  
     clf; n = 1;
[ha, ~] = tight_subplot(3,length(namesP),[0.01 .01],[.08 0],[.05 0]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
   numPoints = 8:size(DataObs_RMSE.(namesP{p}).(glacier),1);

   %Filter data to remove really large RMSE values
   I = DataObs_RMSE.(namesP{p}).(glacier)>1;
   Ftemp = DataObs_RMSE.(namesP{p}).(glacier)(:);
   Ftemp(I(:)) = NaN;
   Ftemp = reshape(Ftemp, size(I));
   
meanWB  = nanmean(Ftemp(8:end,:),2);
stdWB   = nanstd(Ftemp(8:end,:),[],2);

axes(ha(n))
plot(numPoints,meanWB,'LineWidth',0.5,'Color',C(p,:)); hold on
    upper = meanWB + stdWB;
    lower = meanWB - stdWB;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.4,'EdgeColor','none')


    %Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%
        sampledtemp = fullLR.S2.(glacier)(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        estGrid     = diag(sampledtemp);
     RMSEfull = sqrt(mean((estGrid-realGrid.(glacier)(:,1)).^2));

    %Smooth regularily spaced data 
    smoothSize = 7;
    data = DataObs_RMSEeven.(namesP{p}).(glacier)(8:end);
    SS = numPoints(1:size(data,1));
        M = zeros(smoothSize, length(data)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = data(h:end-(smoothSize-h));  end
    T = mean(M);    
    plot(SS((h-1)/2:end-(h-1)/2-1),T,'LineWidth',0.7,'Color','k'); hold on    
     
    %Smooth full data 
    smoothSize = 7;
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = meanWB(h:end-(smoothSize-h));  end
    T = mean(M);    
    
    plot(numPoints((h-1)/2:end-(h-1)/2-1), T,'LineWidth',1.1,'Color',C(p,:)); 
    plot([min(numPoints) max(numPoints)],[RMSEfull,RMSEfull],'--k')

    % 5% of GW mean
    good = find(T-RMSEfull<0.05,1);   
        if ~isempty(good)
        good = good + min(numPoints)-1+(smoothSize-1)/2;
        plot([good good],[0 2],'k:','LineWidth',0.05)
        ConvTable{g,p}  = good;
        end
    %High STD within 25% of mean
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = stdWB(h:end-(smoothSize-h));  end
    TN = mean(M);    
    N10 = find(TN/RMSEfull<0.25,1); 
        if ~isempty(N10) 
        N10 = N10+min(numPoints)-1+(smoothSize-1)/2;
        plot([N10 N10],[0 1.2],'-.k','LineWidth',0.75)
        VarTable{g,p}   = N10;
        end

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
ylim([0 0.6])
xlim([7 45])
n = n+1;

ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);

    if n==2; yInd = 0.002; else yInd = 0.14; end
    if g ==1; sizeG = 0.12; xoff = 0.04; yInd = yInd+0.02; else sizeG = 0.14; xoff = 0.036; end
axes('position',[ax(1)+xoff ax(2)+yInd sizeG sizeG])
pInd = 1:10:length(pUTM.(namesP{p}).(glacier));
if g == 3; fill(options.rig.(glacier)(1:304,1),options.rig.(glacier)(1:304,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
else fill(options.rig.(glacier)(:,1),options.rig.(glacier)(:,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
end
plot(pUTM.(namesP{p}).(glacier)(pInd,1),pUTM.(namesP{p}).(glacier)(pInd,2),'k.','MarkerSize',1.5); 
    axis off; axis equal
end
end 

	%saveFIG_HP('DataObsWB',2,12)
    saveFIG_HP('Pulwicki_Fig5',2,12)
    
    
    
%% Figure 4 - Relative uncertainty
    %clear
%load Patterns.mat 
load Full.mat fullLR
run OPTIONS

    namesP = fieldnames(DataObs_High);  order = [2 3 1 4 5 6]; namesP = namesP(order);
    N1 = {'Centreline',''}; N2 = {'Centreline &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Random',''};
    namesPfull = [N1; N2; N3; N4;N5;N6];


figure(1); clf
SizeToView = 40;
n = 1;

[ha, pos] = tight_subplot(3,length(namesP),[0 0],[0 0],[0 0]);
annStart = [0.061 0.77 .1 .1];

for g = 1:3; 
for p = 1:length(namesP)
    glacier = options.glacier{g};
%Calculate relative uncertainty and plot
    C = fullLR.S2.(glacier);    C = repmat(C,[1 1 nRuns]);
    D = fullSynObs_High.(namesP{p})(SizeToView).(glacier)-C;
    RMSEdist = sqrt(sum(D.^2,3)/nRuns);     
    
        c = nanmean(fullLR.S2.(glacier)(:)); 
        d = nanmean(fullSynObs_High.(namesP{p})(SizeToView).(glacier),2);
        d = nanmean(d,1);   d = d(:)-c;
    RMSEglacier = sqrt(sum(d.^2)/nRuns);
    
axes(ha(n));    
h(n) = imagesc(RMSEdist); hold on
    colormap(cbrewer('seq', 'BuPu', 100,'PCHIP')); set(h(n),'alphadata',~isnan(RMSEdist));
    axis square; axis off;
    caxis([0 0.5])
    
%Plot sampling locations
    E = (UTMinput(SizeToView).(namesP{p}).(glacier)(:,1)-min(options.rig.(glacier)(:,1)))/40;
        minN = min(options.rig.(glacier)(:,2));
        Ng = (options.rig.(glacier)(:,2) - minN)/40; 
    N = (UTMinput(SizeToView).(namesP{p}).(glacier)(:,2)-minN)/40; N = max(Ng)-N;
plot(E,N,'k.','MarkerSize',4.5)

%Add RMSE of glacier-wide WB
    annotation('textbox',[annStart(1)+(p-1)*0.17 annStart(2)-(g-1)*0.3 .1 .1],...
        'String', [num2str(round(RMSEglacier,2),'%1.2f'),' m w.e.'],'EdgeColor','none','FontWeight','bold')

    if g==1; title(namesPfull(p,:)); end

    n = n+1;
end
end

	saveFIG_HP('RelativeUncertainty_SynObs',2,12)
    
    
%% Figure NaN - Elevation Coefficient
%    clear
load Patterns.mat pUTM
load Full.mat fullLR 
load TopoSWE.mat 
run OPTIONS

%%%%%
% Elevation = ElevationReal;    saveName = 'ElevationCoeffReal';    
 Elevation = ElevationSyn;     saveName = 'ElevationCoeffSyn';    
% Elevation = CentreSyn;        saveName = 'CentreCoeffSyn';
% Elevation = SlopeSyn;         saveName = 'SlopeCoeffSyn';
% Elevation = CurveSyn;         saveName = 'CurveCoeffSyn';
% Elevation = WindSyn;          saveName = 'WindCoeffSyn';

    namesP = fieldnames(Elevation);  order = [2 3 1 4 5 7]; namesP = namesP(order);
    N1 = {'Midline',''}; N2 = {'Midline &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Random',''};
    namesPfull = [N1; N2; N3; N4;N5;N6];
    pUTM.Random.G2(:,1) = []; pUTM.Random.G4(:,1) = []; pUTM.Random.G13(:,1) = [];   
    C =[     0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.0588    0.3490    0.1216];

  % Figure  
     clf; n = 1;
     ConvTable   = nan(3,6);  ConvTable = array2table(ConvTable,'VariableNames',namesP);  
     VarTable    = nan(3,6);  VarTable  = array2table(VarTable,'VariableNames',namesP);

[ha, ~] = tight_subplot(3,length(namesP),[0.01 .01],[.08 0],[.05 0]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
   numPoints = 8:size(Elevation.(namesP{p}).(glacier),1);

meanWB   = mean(Elevation.(namesP{p}).(glacier)(8:end,:),2);
stdWB    = std(Elevation.(namesP{p}).(glacier)(8:end,:),[],2);


axes(ha(n));    
plot(numPoints,meanWB,'LineWidth',0.5,'Color',C(p,:)); hold on
    upper = meanWB + stdWB;
    lower = meanWB - stdWB;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.5,'EdgeColor','none'); hold on

%Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%
   
    GWmean = fullLR.S2.coeff{1,g};
    
    smoothSize = 7;
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = meanWB(h:end-(smoothSize-h));  end
    T = mean(M);    
    
    plot(numPoints((h-1)/2:end-(h-1)/2-1), T,'LineWidth',1.1,'Color',C(p,:)); 
    plot([min(numPoints) max(numPoints)],[GWmean,GWmean],'--k')

%     % 5% of GW mean
%     good = find(T-GWmean<0.05,1);   
%         if ~isempty(good)
%         good = good + min(numPoints)-1+(smoothSize-1)/2;
%         plot([good good],[0 2],'k:','LineWidth',0.05)
%         ConvTable{g,p}  = good;
%         end
%     %High STD within 25% of mean
%         M = zeros(smoothSize, length(meanWB)-smoothSize+1);
%     for h = 1:smoothSize;  M(h,:) = stdWB(h:end-(smoothSize-h));  end
%     TN = mean(M);    
%     N10 = find(TN/GWmean<0.25,1); 
%         if ~isempty(N10) 
%         N10 = N10+min(numPoints)-1+(smoothSize-1)/2;
%         plot([N10 N10],[0 1.2],'-.k','LineWidth',0.75)
%         VarTable{g,p}   = N10;
%         end
        

    %titles
    if g==1; title(namesPfull(p,:)); end
    %y labels
    if g ==1 && p==1; ylabel('G4 Elevation Coeff'); set(gca,'XTickLabel',[]);
    elseif g ==2 && p==1; ylabel('G2 Elevation Coeff'); set(gca,'XTickLabel',[]);
    elseif g ==3 && p==1; ylabel('G13 Elevation Coeff'); 
    else set(gca,'YTickLabel',[]);
    end
    if n==15; xlabel('                            Sample size'); end
    if g==1 || g == 2; set(gca,'XTickLabel',[]); end
ylim([-0.1 0.2])
xlim([7 45])
n = n+1;

ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);

    if n==8 || n==9 || n==10 || n==11 || n==12 || n==13 ; yInd = 0.002; else yInd = 0.14; end
    if g ==1; sizeG = 0.12; xoff = 0.04; yInd = yInd+0.02; else sizeG = 0.14; xoff = 0.036; end
axes('position',[ax(1)+xoff ax(2)+yInd sizeG sizeG])
pInd = 1:10:length(pUTM.(namesP{p}).(glacier));
if g == 3; fill(options.rig.(glacier)(1:304,1),options.rig.(glacier)(1:304,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
else fill(options.rig.(glacier)(:,1),options.rig.(glacier)(:,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
end
plot(pUTM.(namesP{p}).(glacier)(pInd,1),pUTM.(namesP{p}).(glacier)(pInd,2),'k.','MarkerSize',1.5); 
    axis off; axis equal
end
end 

    saveFIG_HP(saveName,2,12)