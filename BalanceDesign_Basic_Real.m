%% DATA - WB for all n

% Selecting Data from Pattern

    den = 'S2';
    nRuns = 200;


load TopoSWE.mat
run OPTIONS
clear DataObs_RMSE  

% Remove Accum Area SP
% density_old = [0.348,   0.3327,     0.3487];
% density_new = [0.3422,  0.3436,     0.3692];
%     for g = 1:3;    glacier = options.glacier{g};
%     fullSWE.(den).input.(glacier) = fullSWE.(den).input.(glacier)*density_new(g)/density_old(g);
%     fullSWE.(den).(glacier).swe = fullSWE.(den).(glacier).swe*density_new(g)/density_old(g);
%     end

real_measure    = ObsInCell(fullSWE.(den).input, topo_sampled);

    % Remove dc, aspect and Northness
    for g = 1:3;    glacier = options.glacier{g};
    topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
    
    topo_sampled.(glacier) = rmfield(topo_sampled.(glacier),'centreD');
    topo_sampled.(glacier) = rmfield(topo_sampled.(glacier),'aspect');
    topo_sampled.(glacier) = rmfield(topo_sampled.(glacier),'northness');
    
    topo_sampled_ns.(glacier) = rmfield(topo_sampled_ns.(glacier),'centreD');
    topo_sampled_ns.(glacier) = rmfield(topo_sampled_ns.(glacier),'aspect');
    topo_sampled_ns.(glacier) = rmfield(topo_sampled_ns.(glacier),'northness');
    
    end

load PaperII_AblationArea.mat AblationArea
for g=1:3; glacier = options.glacier{g};
AblationArea.(glacier)(AblationArea.(glacier)==-0.1)=NaN;
AblationArea.(glacier)(~isnan(AblationArea.(glacier)))=1;
end
 
for t = [6,1,3,4,5,100]
if     t == 6; type = 'Circle';           subset = 'pattern';       
elseif t == 1; type = 'Centreline';       subset = 'pattern';     
elseif t == 3; type = 'CentreTransect';   subset = 'pattern';   
elseif t == 4; type = 'Hourglass';        subset = 'pattern';  
elseif t == 5; type = 'HourCircle';       subset = 'pattern';
elseif t == 100; type = 'RandomSafe';     subset = 'random';
end

input.SWE = fullSWE.(den); input.topo_sampled = topo_sampled; 
input.topo_sampled_ns = topo_sampled_ns;

[ subsetSWE_temp, TOPOdata_temp ] = DataSubset( subset, t, input );

[ subsetSWE_temp, TOPOdata_temp ] = ObsInCell( subsetSWE_temp, TOPOdata_temp ); 


 for n = 6:45
     display([type, ' n=',num2str(n)])

for g = 1:3;    glacier = options.glacier{g};

    maxN = length(subsetSWE_temp.(glacier));
    for mc = 1:nRuns

%     nI = randperm(maxN, n);
%     nI = unique(round(linspace(1,length(subsetSWE_temp.(glacier)), n)));
%     WBinput   = subsetSWE_temp.(glacier)(nI,1);
%         ff = fieldnames(TOPOdata_temp.(glacier));
%         TOPOinput = zeros(n,length(ff));
%     for i = 1:length(ff);    fname = ff{i};
%     TOPOinput(:,i)  = TOPOdata_temp.(glacier).(fname)(nI,1); end

%      [WBinput,TOPOinput_temp] = random_uniform(subsetSWE_temp.(glacier), TOPOdata_temp.(glacier), n);
     [WBinput, TOPOinput_temp] = deterministic_locations(subsetSWE_temp.(glacier), TOPOdata_temp.(glacier), n);
        ff = fieldnames(TOPOinput_temp);
        TOPOinput = zeros(n,length(ff));
    for i = 1:length(ff);    fname = ff{i};
    TOPOinput(:,i)  = TOPOinput_temp.(fname); end
    samplingLOC.(type)(n).(glacier) = WBinput(:,2:3);
      
% Linear regression    
        swe	    = WBinput(:,1);
        Xt      = TOPOinput;
        X       = [ones(length(Xt),1), Xt];

        
        % Get coefficients
        coeffs = regress(swe, X);
        MLR(n).(type)(mc).(glacier)   = coeffs;
        
%         display([swe, X])
%         display(coeffs)
        
        %Predict
        sweMLR = repmat(coeffs(1), options.mapsize(g,:));
        mlrCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
            %multiply coeffs and add them
        for m = 1:length(mlrCoeff)
            param               = topoCoeff{m};
            sweT                = topo_full.(glacier).(param)*mlrCoeff(m);
            sweMLR              = sweMLR + sweT;
        end
            %Set min to 0
        sweMLR(sweMLR<0) = 0;
            %Cut accum
%         sweMLR = sweMLR.*AblationArea.(glacier);   
        
%         DataObs.(type).(glacier)(n,mc) = nanmean(sweMLR(:));
         
%         RMSE
        sampledtemp     = sweMLR(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        pred_measure     = diag(sampledtemp);

        real_measure_tmp = real_measure.(glacier)(~isnan(pred_measure));
        pred_measure_tmp = pred_measure(~isnan(pred_measure));
        
        xNaN = ~any(isnan([pred_measure_tmp,real_measure_tmp]),2);
        R = corrcoef(pred_measure_tmp(xNaN), real_measure_tmp(xNaN));
        real_R2.(type).(glacier)(n,mc) = R(2,1)^2;

        realRMSE.(type).(glacier)(n,mc) = sqrt(mean((pred_measure_tmp-real_measure_tmp).^2));
                
    end
        
end

 end
end

% save('PaperII_realdataLR_4var_1000t.mat')
%% Plotting

load PaperII_realdataLR_final density_new density_old
load Patterns.mat pUTM
load PaperII_AblationArea.mat

    namesP = fieldnames(realRMSE);  order = [2 3 1 4 5 6]; namesP = namesP(order);
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
    ela_ind = [1 7; 8 15; 16 23];
    ela_m = [0 0; -0.04*10^5 -0.0088*10^6; 0.04*10^5 0.0086*10^6];


  % Figure  
     clf; n = 1;
[ha, ~] = tight_subplot(3,length(namesP),[0.05 0.01],[.08 0.08],[.07 0.01]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
    
% meanWB  = mean(real_R2.(namesP{p}).(glacier)(numPoints,:),2);
% stdWB   = std(real_R2.(namesP{p}).(glacier)(numPoints,:),[],2);
randuni_meanWB = nanmean(realRMSE.(namesP{p}).(glacier)(numPoints,:),2);
randuni_stdWB  = nanstd(realRMSE.(namesP{p}).(glacier)(numPoints,:),[],2);
det_meanWB = nanmean(det_realRMSE.(namesP{p}).(glacier)(numPoints,:),2);
det_stdWB  = nanstd(det_realRMSE.(namesP{p}).(glacier)(numPoints,:),[],2);
randomWB = random_realRMSE.(namesP{p}).(glacier)(numPoints,1);

% meanR2.(glacier)(numPoints,p) = meanWB;
% stdR2.(glacier)(numPoints,p) = stdWB;
randuni_mean.(glacier)(numPoints,p) = randuni_meanWB;
randuni_std.(glacier)(numPoints,p) = randuni_stdWB;
det_mean.(glacier)(numPoints,p) = det_meanWB;
det_std.(glacier)(numPoints,p) = det_stdWB;
rand_mean.(glacier)(numPoints,p) = randomWB;

axes(ha(n))
plot(numPoints,randuni_meanWB,'LineWidth',0.5,'Color',C(p,:)); hold on
    upper = randuni_meanWB + randuni_stdWB;
    lower = randuni_meanWB - randuni_stdWB;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.4,'EdgeColor','none')
    x2upper = randuni_meanWB + 2*randuni_stdWB;
    x2lower = randuni_meanWB - 2*randuni_stdWB;
fill([numPoints flip(numPoints)],[x2upper',flip(x2lower')],...
     C(p,:),'FaceAlpha',0.2,'EdgeColor','none')

    plot(numPoints,randomWB,'LineWidth',0.5,'Color',[105,105,105]/255); hold on
    plot(numPoints, det_meanWB,':k','LineWidth',1)

    %Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%

    for gg = 1:3; glacier2 = options.glacier{gg};
    alldata.(glacier2) = sweBMS_alldata.(glacier2)*density_new(gg)/density_old(gg); end
        estGrid = SampledCell(alldata);
     best_rmseLR = sqrt(mean((estGrid.(glacier)-realGrid.(glacier)(:,1)).^2));
    plot([min(numPoints) max(numPoints)],[best_rmseLR,best_rmseLR],'--k')
    save('PII_Real_Plot.mat','best_rmseLR','randuni_mean', 'randuni_std', 'det_mean','det_std','rand_mean');

%         R = corrcoef(estGrid.(glacier), realGrid.(glacier)(:,1));
%      bestR2 = R(2,1)^2;
%     plot([min(numPoints) max(numPoints)],[bestR2,bestR2],'--k')
%     save('PII_Real_R2.mat','meanR2','stdR2','bestR2')

%     % 5% of GW mean
% %     good = find(T-bestRMSE<0.05,1);  
%     good = find(meanWB-bestRMSE<0.05,1);
% %       good = find(2*stdWB<0.05,1);
%         if ~isempty(good)
%         good = good + min(numPoints);%-1+(smoothSize-1)/2;
% %         plot([good good],[0 2],'k-.','LineWidth',0.05)
%         ConvTable{g,p}  = good;
%         end

    %titles
    if g==1; title(namesPfull(p,:)); end
    %y labels
    if g ==1 && p==1; ylabel({'G4 RMSE', '(m w.e.)'}); set(gca,'XTickLabel',[]);
    elseif g ==2 && p==1; ylabel({'G2 RMSE', '(m w.e.)'}); set(gca,'XTickLabel',[]);
    elseif g ==3 && p==1; ylabel({'G13 RMSE', '(m w.e.)'}); 
    else set(gca,'YTickLabel',[]);
    end
    if n==15; xlabel('                            Sample size (N)'); end
    if g==1 || g == 2; set(gca,'XTickLabel',[]); end
ylim([0 0.4])
xlim([min(numPoints) max(numPoints)])
fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)


ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);


    if n==0; yInd = 0.002; else; yInd = 0.13; end
%     yInd = 0.14;
    if g ==1; sizeG = 0.12; xoff = 0.035; yInd = yInd+0.02; else sizeG = 0.14; xoff = 0.03; end
axes('position',[ax(1)+xoff ax(2)+yInd sizeG sizeG])
pInd = 1:length(pUTM.(namesP{p}).(glacier));
if g == 3; fill(options.rig.(glacier)(1:304,1),options.rig.(glacier)(1:304,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
else; fill(options.rig.(glacier)(:,1),options.rig.(glacier)(:,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
end
plot(pUTM.(namesP{p}).(glacier)(pInd,1),pUTM.(namesP{p}).(glacier)(pInd,2),'k.','MarkerSize',1.5); 
plot(options.rig.ELA(ela_ind(g,1):ela_ind(g,2),1)+ela_m(g,1), options.rig.ELA(ela_ind(g,1):ela_ind(g,2),2)+ela_m(g,2),'-k');
    axis off; axis equal
n = n+1;

end
% display(bestRMSE)
end 

%     saveFIG_HP('PII_AA_RealData_wRandom',2,12)
    
%% nc and nv calculation

% load PII_FastRuns.mat
load PaperII_AblationArea.mat sweBMS_alldata
load Full.mat fullLR fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS

numPoints   = 6:45;
namesP      = fieldnames(realRMSE);
nc_table    = nan(3,6);     nc_table = array2table(nc_table,'VariableNames',namesP);  
nv_table    = nan(3,6);     nv_table = array2table(nv_table,'VariableNames',namesP);  
realGrid    = ObsInCell(fullSWE.S2.input, topo_sampled);

for g = 1:3; glacier = options.glacier{g};
    estGrid     = SampledCell(sweBMS_alldata);
    bestRMSE    = sqrt(mean((estGrid.(glacier)-realGrid.(glacier)(:,1)).^2));
for p = 1:length(namesP)
    
meanWB  = mean(realRMSE.(namesP{p}).(glacier)(numPoints,:),2);
stdWB   = std(realRMSE.(namesP{p}).(glacier)(numPoints,:),[],2);

    %Smooth mean 
    smoothSize = 3;
        M = zeros(smoothSize, length(meanWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = meanWB(h:end-(smoothSize-h));  end
    meanSmooth = mean(M); 
    %Smooth std  
        M = zeros(smoothSize, length(stdWB)-smoothSize+1);
    for h = 1:smoothSize;  M(h,:) = stdWB(h:end-(smoothSize-h));  end
    stdSmooth = mean(M);    



nc = find(meanSmooth-bestRMSE<bestRMSE*0.1,1);
if isempty(nc); nc = nan; 
else; nc = nc + min(numPoints); end
nc_table{g,p}   = nc;

nv = find(stdSmooth+meanSmooth-bestRMSE<bestRMSE*0.25,1);
if isempty(nv); nv=nan; 
else; nv = nv + min(numPoints); end
nv_table{g,p}   = nv;
    
end 
end