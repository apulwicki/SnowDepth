%% DATA - WB for all n

% Selecting Data from Pattern

    den = 'S2';
    nRuns = 100;


load TopoSWE.mat
run OPTIONS
clear DataObs_RMSE  

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

    nI = randperm(maxN, n);
    WBinput   = subsetSWE_temp.(glacier)(nI,1);
        ff = fieldnames(TOPOdata_temp.(glacier));
        TOPOinput = zeros(n,length(ff));
    for i = 1:length(ff);    fname = ff{i};
    TOPOinput(:,i)  = TOPOdata_temp.(glacier).(fname)(nI,:); end
      
% Linear regression    
        swe	    = WBinput;
        Xt      = TOPOinput;
        X       = [ones(length(Xt),1), Xt];

        % Get coefficients
        coeffs = regress(swe, X);
        MLR(n).(type)(mc).(glacier)   = coeffs;

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
        sweMLR = sweMLR.*AblationArea.(glacier);   
        
%         DataObs.(type).(glacier)(n,mc) = nanmean(sweMLR(:));
         
%         RMSE
        sampledtemp     = sweMLR(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        pred_measure     = diag(sampledtemp);

        real_measure_tmp = real_measure.(glacier)(~isnan(pred_measure));
        pred_measure_tmp = pred_measure(~isnan(pred_measure));

        realRMSE.(type).(glacier)(n,mc) = sqrt(mean((pred_measure_tmp-real_measure_tmp).^2));
                
    end
        
end

 end
end

%% Plotting

load PII_FastRuns.mat
load Full.mat fullLR fullSWE
load TopoSWE.mat topo_sampled
run OPTIONS
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
[ha, ~] = tight_subplot(3,length(namesP),[0.01 .01],[.08 0],[.05 0]);
for g = 1:3; glacier = options.glacier{g};
for p = 1:length(namesP)
%    numPoints = 8:size(DataObs_RMSE.(namesP{p}).(glacier),1);
    
meanWB  = mean(realRMSE.(namesP{p}).(glacier)(numPoints,:),2);
stdWB   = std(realRMSE.(namesP{p}).(glacier)(numPoints,:),[],2);


axes(ha(n))
plot(numPoints,meanWB,'LineWidth',0.5,'Color',C(p,:)); hold on
    upper = meanWB + stdWB;
    lower = meanWB - stdWB;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.4,'EdgeColor','none')


    %Best sample size Fitting function %%%%%%%%%%%%%%%%%%%%%%%%%%

        estGrid = SampledCell(sweBMS_alldata);
     RMSEfull = sqrt(mean((estGrid.(glacier)-realGrid.(glacier)(:,1)).^2));

 
    plot([min(numPoints) max(numPoints)],[RMSEfull,RMSEfull],'--k')

    % 5% of GW mean
%     good = find(T-RMSEfull<0.05,1);  
    good = find(meanWB-RMSEfull<0.05,1);
%       good = find(2*stdWB<0.05,1);
        if ~isempty(good)
        good = good + min(numPoints);%-1+(smoothSize-1)/2;
        plot([good good],[0 2],'k-.','LineWidth',0.05)
        ConvTable{g,p}  = good;
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
ylim([0 0.5])
xlim([min(numPoints) 25])
n = n+1;

ax = get(gca,'Position');
set(gca,'Position', [ax(1) ax(2) ax(3) 0.28]);

%     if n==2; yInd = 0.002; else; yInd = 0.14; end
    yInd = 0.14;
    if g ==1; sizeG = 0.12; xoff = 0.04; yInd = yInd+0.02; else sizeG = 0.14; xoff = 0.036; end
axes('position',[ax(1)+xoff ax(2)+yInd sizeG sizeG])
pInd = 1:10:length(pUTM.(namesP{p}).(glacier));
if g == 3; fill(options.rig.(glacier)(1:304,1),options.rig.(glacier)(1:304,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
else; fill(options.rig.(glacier)(:,1),options.rig.(glacier)(:,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
end
plot(pUTM.(namesP{p}).(glacier)(pInd,1),pUTM.(namesP{p}).(glacier)(pInd,2),'k.','MarkerSize',1.5); 
plot(options.rig.ELA(ela_ind(g,1):ela_ind(g,2),1)+ela_m(g,1), options.rig.ELA(ela_ind(g,1):ela_ind(g,2),2)+ela_m(g,2),'-.k');
    axis off; axis equal
end
end 

%     saveFIG_HP('PII_AA_RealData',2,12)
    