%% Predict bw and calculate RMSE
run OPTIONS.m
load TestSFU.mat BMS_HourCircle
load TopoSWE.mat fullSWE topo_sampled topo_full
load PaperII_AblationArea.mat AblationArea


% ind = 6:3:45;

    % Remove dc, aspect and Northness
for g = 1:3;    glacier = options.glacier{g};
    topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
end

% Get coeffs
ind = 6;

for num_model = 1:size(BMS_HourCircle.G4,2)
    
for g = 1:3;    glacier = options.glacier{g};
coeffs.(glacier) = [BMS_HourCircle.(glacier)(ind,num_model,1), BMS_HourCircle.(glacier)(ind,num_model,2),...
                    BMS_HourCircle.(glacier)(ind,num_model,3), BMS_HourCircle.(glacier)(ind,num_model,4),...
                    BMS_HourCircle.(glacier)(ind,num_model,5)];
end

real_measure = ObsInCell(fullSWE.S2.input, topo_sampled);

% patterns  = fieldnames(BMS);
% 
% for p = 1:6; pattern_name = patterns{p};
%     est_Bw.(pattern_name) = zeros(length(BMS),3);
    
% for r = 5:length(BMS)
    
for g = 1:3;    glacier = options.glacier{g};
    sweBMS      = repmat(coeffs.(glacier)(5), options.mapsize(g,:));
    mlrCoeff    = coeffs.(glacier)(1:4);    topoCoeff = fieldnames(topo_full.G4);
        %multiply coeffs and add them
    for m = 1:length(mlrCoeff)
        param     = topoCoeff{m};
        sweT      = topo_full.(glacier).(param)*mlrCoeff(m);
        sweBMS    = sweBMS + sweT;
    end
        %Set min to 0
    sweBMS(sweBMS<0) = 0;

    % Nan all spots where accum area
    sweBMS(AblationArea.(glacier)==-0.1) = NaN;
    pred_bw.(glacier) = sweBMS;
end

    % Calculate RMSE
        pred_measure = SampledCell(pred_bw);
        
        for gg = 1:3;        glacier = char(options.glacier(gg));
        real_measure_tmp = real_measure.(glacier)(~isnan(pred_measure.(glacier)));
        pred_measure_tmp = pred_measure.(glacier)(~isnan(pred_measure.(glacier)));
        
        realRMSE.(namesP{p}).(glacier)(ind,num_model) = sqrt(mean((pred_measure_tmp-real_measure_tmp).^2));
        end

end

% %% Mean and std of Real RMSE
% for gg = 1:3;        glacier = char(options.glacier(gg));
%    mean_realRMSE.(glacier)(ind)    = mean(realRMSE.(glacier)(ind,:));
%    std_realRMSE.(glacier)(ind)     = std(realRMSE.(glacier)(ind,:));
% end

%% Figure 5 - Real data
load Full.mat fullLR fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS

    namesP = fieldnames(realRMSE);  order = [2 3 1 4 5 6]; namesP = namesP(order);
    N1 = {'Midline',''}; N2 = {'Mid &','Transverse'}; N3 = {'Circle',''};
    N4 = {'Hourglass',''}; N5 = {'Hourglass','& Circle'}; N6 = {'Random',''}; 
    namesPfull = [N1; N2; N3; N4;N5;N6];
%     pUTM.Random.G2(:,1) = []; pUTM.Random.G4(:,1) = []; pUTM.Random.G13(:,1) = [];   
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
%    numPoints = 8:size(DataObs_RMSE.(namesP{p}).(glacier),1);
    numPoints = 6:3:45;
    
   %Filter data to remove really large RMSE values
   I = realRMSE.(namesP{p}).(glacier)>1;
   Ftemp = realRMSE.(namesP{p}).(glacier)(:);
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
    data = realRMSE.(namesP{p}).(glacier)(8:end);
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
%     saveFIG_HP('Pulwicki_Fig5',2,12)
    
