%% Figure 2 - Synthetic sampling designs 

    load Patterns.mat pUTM 
    load Full.mat fullLR
    run OPTIONS
NumSubPoints = 15;
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
PlotTopoParameter_DD(fullLR.S2, 'WB (m w.e.)', patternFULL, patternSUB, 'black')
	saveFIG_HP(['SampleDesign_',P{t}],1,8.6)
end

%% Figure 3 - Synthetic data WB 
    clear
load Patterns.mat 
load Full.mat fullLR fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS

    namesP = fieldnames(SynObs_High);  order = [2 3 1 4 5 7]; namesP = namesP(order);
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
[ha, ~] = tight_subplot(3,length(namesP),[0.01 .01],[.08 0],[.05 0]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
   numPoints = 8:size(SynObs_High.(namesP{p}).(glacier),1);

    M1   = mean(SynObs_High.(namesP{p}).(glacier)(8:end,:),2);
    M2   = mean(SynObs_Low.(namesP{p}).(glacier)(8:end,:),2);   
meanWB   = mean([M1 M2],2);
stdWB_H  = std(SynObs_High.(namesP{p}).(glacier)(8:end,:),[],2);
stdWB_L  = std(SynObs_Low.(namesP{p}).(glacier)(8:end,:),[],2);


axes(ha(n));    
plot(numPoints,meanWB,'LineWidth',1.1,'Color',C(p,:)); hold on
    upper = meanWB + stdWB_H;
    lower = meanWB - stdWB_H;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.2,'EdgeColor','none')
     upper = meanWB + stdWB_L;
    lower = meanWB - stdWB_L;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.4,'EdgeColor','none')

%Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%
    
    GWmean = nanmean(fullLR.S2.(glacier)(:));
    F = fit(numPoints',meanWB,'exp2');
    T = F(numPoints);
    
    plot(numPoints, T,'k'); hold on
    plot([min(numPoints) max(numPoints)],[GWmean,GWmean],'--k')

    % 5% of GW mean
    good = find(T-GWmean<0.05,1);
        if ~isempty(good)
        plot([good+min(numPoints) good+min(numPoints)],[0 2],'k:','LineWidth',0.05)
        end
    % 1% of GW mean
%     good = find(T-GWmean<0.01*GWmean,1);
%         if ~isempty(good)
%         plot([good+min(numPoints) good+min(numPoints)],[0 2],'k-.','LineWidth',0.25)
%         end
    %High STD within 10% of mean
    FN  = fit(numPoints',stdWB_H,'exp2');
    TN  = FN(numPoints);
    N10 = find(TN<0.15,1); N10 = N10+min(numPoints)-1;
    %N10 = find(stdWB_H<0.1,1);  N10 = N10+min(numPoints)-1;
        if ~isempty(N10) 
        plot([N10 N10],[0 1.2],'-.k','LineWidth',0.75)
        end

    %titles
    if g==1; title(namesPfull(p,:)); end
    %y labels
    if g ==1 && p==1; ylabel('G4 WB (m w.e.)'); set(gca,'XTickLabel',[]);
    elseif g ==2 && p==1; ylabel('G2 WB (m w.e.)'); set(gca,'XTickLabel',[]);
    elseif g ==3 && p==1; ylabel('G13 WB (m w.e.)'); 
    else set(gca,'YTickLabel',[]);
    end
    if n==15; xlabel('                    Sample size'); end
    if g==1 || g == 2; set(gca,'XTickLabel',[]); end
ylim([0.2 1.1])
xlim([7 45])
n = n+1;

ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);

    if n==1 || n==8 || n == 14; yInd = 0.002; else yInd = 0.14; end
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

	saveFIG_HP('SyntheticObsWB',2,12)

%% Figure 4 - Real data
    clear
load Patterns.mat 
load Full.mat fullLR fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS

    namesP = fieldnames(DataObs_HighRMSE);  order = [2 3 1 4 5 6]; namesP = namesP(order);
    N1 = {'Midline',''}; N2 = {'Mid &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Safe',''}; 
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
[ha, ~] = tight_subplot(3,length(namesP),[0.01 .01],[.08 0],[.05 0]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
   numPoints = 8:size(DataObs_RMSE.(namesP{p}).(glacier),1);

meanWB   = mean(DataObs_RMSE.(namesP{p}).(glacier)(8:end,:),2);
%stdWB_H  = std(DataObs_HighRMSE.(namesP{p}).(glacier)(8:L+7,:),[],2);
%stdWB_L  = std(DataObs_LowRMSE.(namesP{p}).(glacier)(8:L+7,:),[],2);

axes(ha(n))
plot(numPoints,meanWB,'LineWidth',1.1,'Color',C(p,:)); hold on
%     upper = meanWB + stdWB_H;
%     lower = meanWB - stdWB_H;
% fill([numPoints flip(numPoints)],[upper',flip(lower')],...
%      C(p,:),'FaceAlpha',0.2,'EdgeColor','none')
%     upper = meanWB + stdWB_L;
%     lower = meanWB - stdWB_L;
% fill([numPoints flip(numPoints)],[upper',flip(lower')],...
%      C(p,:),'FaceAlpha',0.4,'EdgeColor','none')


    %Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%
    F = fit(numPoints',meanWB,'exp2');
    T = F(numPoints);
        sampledtemp = fullLR.S2.(glacier)(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        estGrid     = diag(sampledtemp);
        realGrid    = ObsInCell(fullSWE.S2.input, topo_sampled);
    RMSEfull = sqrt(mean((estGrid-realGrid.(glacier)(:,1)).^2));

    plot(numPoints, T,'k'); hold on
    plot([min(numPoints) max(numPoints)],[RMSEfull,RMSEfull],'--k')

    % difference of <0.05
    good = find(T-RMSEfull<0.05,1);
        if ~isempty(good)
        plot([good+min(numPoints) good+min(numPoints)],[0 2],'k:','LineWidth',0.05)
        end
    % 50% of GW mean
%     good = find(T-RMSEfull<0.50*RMSEfull,1);
%         if ~isempty(good)
%         plot([good+min(numPoints) good+min(numPoints)],[0 2],'k-.','LineWidth',0.25)
%         end
    %High STD <0.01
%     FN  = fit(numPoints',stdWB_H,'exp2');
%     TN  = FN(numPoints);
%     N10 = find(TN<0.01,1); N10 = N10+min(numPoints)-1;
%     %N10 = find(stdWB_H<0.1,1);  N10 = N10+min(numPoints)-1;
%         if ~isempty(N10) 
%         plot([N10 N10],[0 1.2],'-.k','LineWidth',0.5)
%         end


    %titles
    if g==1; title(namesPfull(p,:)); end
    %y labels
    if g ==1 && p==1; ylabel('G4 RMSE (m w.e.)'); set(gca,'XTickLabel',[]);
    elseif g ==2 && p==1; ylabel('G2 RMSE (m w.e.)'); set(gca,'XTickLabel',[]);
    elseif g ==3 && p==1; ylabel('G13 RMSE (m w.e.)'); 
    else set(gca,'YTickLabel',[]);
    end
    if n==15; xlabel('                    Sample size'); end
    if g==1 || g == 2; set(gca,'XTickLabel',[]); end
ylim([0 0.6])
xlim([7 45])
n = n+1;

ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);

    if n==2 || n == 14; yInd = 0.002; else yInd = 0.14; end
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

	saveFIG_HP('DataObsWB',2,12)

    
%% Figure 5 - Relative uncertainty
    %clear
load Patterns.mat 
load Full.mat fullLR
run OPTIONS

    namesP = fieldnames(DataObs_High);  order = [2 3 1 4 5 6]; namesP = namesP(order);
    N1 = {'Midline',''}; N2 = {'Mid &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Safe',''};
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
    caxis([0 1])
    
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
