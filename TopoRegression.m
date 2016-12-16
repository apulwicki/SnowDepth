%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

load TopoParams.mat
%run Import_Topo.m

%% MLR - Topo Regression

% %remove aspect
% glacier = {'G4','G2','G13'};
% for i = 1:3
%    name     = char(glacier(i));
%    topo_sampled.(name) = rmfield(topo_sampled.(name),'aspect');
% end

for t = 2:9
run OPTIONS.m
options.DensitySWE  = t;
options.ZZ          = 2; %exclude zigzags
run MAIN

    for i = 1:3
        y       = SWE(i).swe;
        glacier = char(options.glacier(i)); 
            display(['option = ',num2str(t), ', glacier = ',glacier]);
        X       = topo_sampled.(glacier);

        [MLR(t).(glacier), residualsMLR(t).(glacier)] = MLRcalval(y, X);
        MLR(t).(glacier).Properties.VariableNames = {['MLRCoefficient', num2str(t)],['MLRPercentVarExplaned', num2str(t)]};
    end
end
        clear best i name X y t glacier
        
%% Export all values
G4_mlrDensity = []; G2_mlrDensity = []; G13_mlrDensity = [];
for i = 2:9
G4_mlrDensity  = [G4_mlrDensity, MLR(i).G4(:,3)];
G2_mlrDensity  = [G2_mlrDensity, MLR(i).G2(:,3)];
G13_mlrDensity = [G13_mlrDensity, MLR(i).G13(:,3)];
end
%     writetable(G4_mlrDensity,'/home/glaciology1/Documents/Data/GlacierTopos/G4_mlrDensity.xlsx')
%     writetable(G2_mlrDensity,'/home/glaciology1/Documents/Data/GlacierTopos/G2_mlrDensity.xlsx')
%     writetable(G13_mlrDensity,'/home/glaciology1/Documents/Data/GlacierTopos/G13_mlrDensity.xlsx')
%     
    writetable(G4_mlrDensity,'/Users/Alexandra/Documents/SFU/Data/G4_mlrDensity.xlsx')
    writetable(G2_mlrDensity,'/Users/Alexandra/Documents/SFU/Data/G2_mlrDensity.xlsx')
    writetable(G13_mlrDensity,'/Users/Alexandra/Documents/SFU/Data/G13_mlrDensity.xlsx')
    
    
%Check normality of reisduals
%plotResiduals(lm.G4) %many options for plotting

clear A box centreline d distance f* G* h i min_dist NoZZ param* RGB X1 Y1

%% ANOVA between topographic params

for i = 1:3
    y  = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    topo    = [aspect.(name); northness.(name); profileCurve.(name); ...
                tangentCurve.(name); slope.(name); elevation.(name); Sx.(name)];
    group   = cellstr([repmat('aspect',length(aspect.(name)),1); repmat('northn',length(northness.(name)),1);...
                repmat('profCu',length(profileCurve.(name)),1); repmat('tangCu',length(tangentCurve.(name)),1);...
                repmat('slopee',length(slope.(name)),1); repmat('elevat',length(elevation.(name)),1); repmat('Sxxxxx',length(Sx.(name)),1)]);
    [p,tbl,stats] = anova1(topo,group);
        
            
end    
    
%% Correlation between topographic parameters 

%Sampled
for i = 1:3
X = [];
glacier = char(options.glacier(i));
params = fieldnames(topo_sampled.(glacier));
    for p = 1:length(params)
        PP = char(params(p));
        X = [X, topo_sampled.(glacier).(PP)];
    end
    [pearson.(glacier) Ppearson.(glacier)] = corr(X); 
    pearson.(glacier)   = round(triu(pearson.(glacier)),2); 
    Ppearson.(glacier)  = round(triu(Ppearson.(glacier)),2);
%     matrix2latex(pearson.(glacier),'/home/glaciology1/Documents/MastersDocuments/Methods/temp.txt',...
%         'rowLabels',options.topoVars_xunit,'columnLabels',options.topoVars_xunit, 'alignment','c')
end
    
%Full topo
for i = 1:3
X = [];
glacier = char(options.glacier(i));
params = fieldnames(topo_full.(glacier));
    for p = 1:length(params)
        PP = char(params(p));
        X = [X, topo_full.(glacier).(PP)(:)];
    end
    [pearson.(glacier) Ppearson.(glacier)] = corr(X); 
    pearson.(glacier)   = round(triu(pearson.(glacier)),2); 
    Ppearson.(glacier)  = round(triu(Ppearson.(glacier)),2);
end
    
%% BMS

for t = 2:9
    run OPTIONS
    options.DensitySWE  = t;
    options.ZZ          = 2; %exclude zigzags
    run MAIN
    display(['Option ', num2str(t)]);
for g = 1:3
    glacier = char(options.glacier(g));
    
    cd BMS
    [BMSinit, BMSres] = BMS_R(SWE, topo_sampled);
    cd ..
BMS(t).(glacier)     = BMSinit.(glacier);   
BMS(t).(glacier).Properties.VariableNames = {['BMSCoefficient', num2str(t)],['BMSPercentVarExplaned', num2str(t)]};

residualsBMS(t).(glacier) = BMSres.(glacier);
end
end
    clear BMSinit BMSres t g glacier


%% Predicting

for t = 2:9
    %check coeff order - MLR
    mlrCoeff = MLR(t).G4.Properties.RowNames(1:end-2);   topoCoeff = fieldnames(topo_full_ns.G4);
    if ~isequal(mlrCoeff, topoCoeff)
        disp('Different order of coefficients between MLR and topo'); return; end
    %check coeff order - BMS
    bmaCoeff = BMS(t).G4.Properties.RowNames(1:end-2);   topoCoeff = fieldnames(topo_full_ns.G4);
    if ~isequal(bmaCoeff, topoCoeff)
        disp('Different order of coefficients between BMS and topo'); return; end
    
    for g = 1:3
    glacier = char(options.glacier(g));
        %MLR
         %Intercept
        sweMLR(t).(glacier) = repmat(MLR(t).(glacier){end-1,1}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        for n = 1:length(mlrCoeff)
            param               = char(mlrCoeff(n));
            sweT                = topo_full.(glacier).(param)*MLR(t).(glacier){n,1};
            sweMLR(t).(glacier) = sweMLR(t).(glacier) + sweT;
        end
        
        %BMS
         %Intercept
        sweBMS(t).(glacier) = repmat(BMS(t).(glacier){end-1,1}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        for n = 1:length(bmaCoeff)
            param               = char(bmaCoeff(n));
            sweT                = topo_full.(glacier).(param)*BMS(t).(glacier){n,1};
            sweBMS(t).(glacier) = sweBMS(t).(glacier) + sweT;
        end
    end
end

    clear param sweT *Coeff glacier g n t

%% Return min and max regression coeffs

%Rearrange to compare density options
for g = 1:3
    glacier = char(options.glacier(g));
    boxVarBMS.(glacier) = [];
    boxVarMLR.(glacier) = [];
for i = 2:9
        BMS(i).(glacier).Properties.VariableNames(1,1) = {['Option', num2str(i-1)]};
    boxVarBMS.(glacier)  = [boxVarBMS.(glacier),  BMS(i).(glacier)(1:end-2,2)];%,  MLR(i).G4(1:end-1,2)];
    
        MLR(i).(glacier).Properties.VariableNames(1,1) = {['Option', num2str(i-1)]};
    boxVarMLR.(glacier)  = [boxVarMLR.(glacier),    MLR(i).(glacier)(1:end-2,2)];
end
        boxVarBMS.(glacier).Properties.RowNames = strcat(boxVarBMS.(glacier).Properties.RowNames,glacier);
        boxVarMLR.(glacier).Properties.RowNames = strcat(boxVarMLR.(glacier).Properties.RowNames,glacier);
end

%Save to table
    writetable([boxVarMLR.G4;boxVarMLR.G2;boxVarMLR.G13] ,'/home/glaciology1/Downloads/MLR')
    writetable([boxVarBMS.G4;boxVarBMS.G2;boxVarBMS.G13] ,'/home/glaciology1/Downloads/BMS')


%boxplot(boxVarMLR.G4{1:end-1,:}','labels',options.topoVars)

for g = 1:3
    glacier = char(options.glacier(g));
RegressC.(glacier) = table(mean(box.(glacier){:,:},2), min(box.(glacier){:,:},[],2), max(box.(glacier){:,:},[],2),...
                'VariableNames',{'Mean', 'Min','Max'}, 'RowNames',BMS(2).G4.Properties.RowNames(1:end-1));
end
    
    clear g i box glacier






