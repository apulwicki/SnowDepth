%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

load('TopoSWE.mat')
%% MLR - Topo Regression

% %remove aspect
% glacier = {'G4','G2','G13'};
% for i = 1:3
%    name     = char(glacier(i));
%    topo_sampled.(name) = rmfield(topo_sampled.(name),'aspect');
% end

for t = 2:9
    for i = 1:3
        glacier = char(options.glacier(i)); 
        swe       = sweOPT(t).(glacier)(:,1);
            display(['option = ',num2str(t), ', glacier = ',glacier]);
        X       = topo_sampled.(glacier);

        [MLR(t).(glacier), residualsMLR(t).(glacier)] = MLRcalval(swe, X);
        MLR(t).(glacier).Properties.VariableNames = {['MLRCoefficient_', num2str(t)],...
                                                     ['MLRsemiR2_', num2str(t)],...
                                                     ['MLRunivarR2_', num2str(t)]};
    end
end

for g = 1:3
    glacier = char(options.glacier(g)); 
        for i = 2:9
        stackMLR(:,:,i-1) = MLR(i).(glacier){:,:}; end
        meanMLR = mean(stackMLR, 3);

MLR(10).(glacier) = table(meanMLR(:,1),meanMLR(:,2),meanMLR(:,3),...
                    'RowNames',MLR(9).G4.Properties.RowNames,...
                    'VariableNames',{'MLRmeanCoeff','MLRmeanSemiR2','MLRmeanUniR2'});
end

display('Done');

        clear best i name X y t glacier stackMLR meanMLR g
        
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
    swe  = SWE(i).swe;
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
glacier = char(options.glacier(i));    
X = struct2table(topo_sampled_ns.(glacier));
    
    [pearson.(glacier) Ppearson.(glacier)] = corr(X{:,:}); 
    pearson.(glacier)   = round(triu(pearson.(glacier)),2); 
    Ppearson.(glacier)  = round(triu(Ppearson.(glacier)),2);
    
%     pearsonT.(glacier)  = table(pearson.(glacier)(:,1),pearson.(glacier)(:,2),...
%                                 pearson.(glacier)(:,3),pearson.(glacier)(:,4),...
%                                 pearson.(glacier)(:,5),pearson.(glacier)(:,6),...
%                                 pearson.(glacier)(:,7),pearson.(glacier)(:,8),...
%                             'RowNames', options.topoVars);
%     writetable(pearsonT.(glacier),['/home/glaciology1/Downloads/corr',glacier])
end
    
%Full topo
for i = 1:3
glacier = char(options.glacier(i));
X = struct2table(topo_sampled_ns.(glacier));

    [pearson.(glacier) Ppearson.(glacier)] = corr(X); 
    pearson.(glacier)   = round(triu(pearson.(glacier)),2); 
    Ppearson.(glacier)  = round(triu(Ppearson.(glacier)),2);
end
    
%% BMA

for t = 2:9
    swe       = sweOPT(t);
        display(['option = ',num2str(t)]);

    cd BMS
    [BMSinit, BMSres] = BMS_R(swe, topo_sampled);
    cd ..
    
for g = 1:3
    glacier = char(options.glacier(g));
BMS(t).(glacier) = BMSinit.(glacier);   
BMS(t).(glacier).Properties.VariableNames = {['BMSCoefficient_', num2str(t)],...
                                             ['BMSsemiR2_', num2str(t)],...
                                             ['BMSunivarR2_', num2str(t)]};
residualsBMS(t).(glacier) = BMSres.(glacier);
end
end

for g = 1:3
    glacier = char(options.glacier(g)); 
        for i = 2:9
        stackMLR(:,:,i-1) = BMS(i).(glacier){:,:}; end
        meanMLR = mean(stackMLR, 3);

BMS(10).(glacier) = table(meanMLR(:,1),meanMLR(:,2),meanMLR(:,3),...
                    'RowNames',BMS(9).G4.Properties.RowNames,...
                    'VariableNames',{'BMSmeanCoeff','BMSmeanSemiR2','BMSmeanUniR2'});
end

display('Done');

        clear best i name X y t glacier stackMLR meanMLR g

    clear BMSinit BMSres t g glacier swe


%% Predicting

for t = 2:9
    %check coeff order - MLR
    mlrCoeff = MLR(t).G4.Properties.RowNames(1:end-3);   topoCoeff = fieldnames(topo_full_ns.G4);
    if ~isequal(mlrCoeff, topoCoeff)
        disp('Different order of coefficients between MLR and topo'); return; end
    %check coeff order - BMS
    bmaCoeff = BMS(t).G4.Properties.RowNames(1:end-3);   topoCoeff = fieldnames(topo_full_ns.G4);
    if ~isequal(bmaCoeff, topoCoeff)
        disp('Different order of coefficients between BMS and topo'); return; end
    
    for g = 1:3
    glacier = char(options.glacier(g));
        %MLR
         %Intercept
        sweMLR(t).(glacier) = repmat(MLR(t).(glacier){end-2,1}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        for n = 1:length(mlrCoeff)
            param               = char(mlrCoeff(n));
            sweT                = topo_full.(glacier).(param)*MLR(t).(glacier){n,1};
            sweMLR(t).(glacier) = sweMLR(t).(glacier) + sweT;
        end
         %Set min to 0
        sweMLR(t).(glacier)(sweMLR(t).(glacier)<0) = 0;
        
        %BMS
         %Intercept
        sweBMS(t).(glacier) = repmat(BMS(t).(glacier){end-2,1}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        for n = 1:length(bmaCoeff)
            param               = char(bmaCoeff(n));
            sweT                = topo_full.(glacier).(param)*BMS(t).(glacier){n,1};
            sweBMS(t).(glacier) = sweBMS(t).(glacier) + sweT;
        end
         %Set min to 0
        sweBMS(t).(glacier)(sweBMS(t).(glacier)<0) = 0;
    end
end

    clear param sweT *Coeff glacier g n t

for g = 1:3
    glacier = char(options.glacier(g));
    
    for i = 2:9
    stackSWE.MLR.(glacier)(:,:,i-1)   = sweMLR(i).(glacier);
    stackSWE.BMS.(glacier)(:,:,i-1)   = sweBMS(i).(glacier);
    end

end    
    clear g i glacier
    
%% Return min and max regression coeffs

%Rearrange to compare density options
for g = 1:3
    glacier = char(options.glacier(g));
    boxBMS.(glacier) = [];
    boxMLR.(glacier) = [];
    boxALL.(glacier)    = [];
for i = 2:9
        BMS(i).(glacier).Properties.VariableNames(1,1) = {['BMSOption', num2str(i-1)]};
    boxBMS.(glacier)  = [boxBMS.(glacier),  BMS(i).(glacier)(1:end-1,1)];%,  MLR(i).G4(1:end-1,2)];
    
        MLR(i).(glacier).Properties.VariableNames(1,1) = {['MLROption', num2str(i-1)]};
    boxMLR.(glacier)  = [boxMLR.(glacier),    MLR(i).(glacier)(1:end-1,1)];
        
    boxALL.(glacier)    = [boxALL.(glacier),  BMS(i).(glacier)(1:end-1,1),  MLR(i).(glacier)(1:end-1,1)];

end
        boxBMS.(glacier).Properties.RowNames = strcat(boxBMS.(glacier).Properties.RowNames,glacier);
        boxMLR.(glacier).Properties.RowNames = strcat(boxMLR.(glacier).Properties.RowNames,glacier);
end

%Save to table
%     writetable([boxMLR.G4;boxMLR.G2;boxMLR.G13] ,'/home/glaciology1/Downloads/MLR')
%     writetable([boxBMS.G4;boxBMS.G2;boxBMS.G13] ,'/home/glaciology1/Downloads/BMS')

for g = 1:3
    glacier = char(options.glacier(g));
RegressC.(glacier) = table(mean(boxALL.(glacier){:,:},2), min(boxALL.(glacier){:,:},[],2), max(boxALL.(glacier){:,:},[],2),...
                'VariableNames',{'Mean', 'Min','Max'}, ...
                'RowNames',BMS(2).G4.Properties.RowNames(1:end-1));
end
%     writetable([RegressC.G4;RegressC.G2;RegressC.G13] ,'/home/glaciology1/Downloads/Coeffs')
% 'RowNames',strcat(BMS(2).G4.Properties.RowNames(1:end-1),glacier))
    clear g i glacier

%Predict swe with new coeffs    
for g = 1:3
    glacier = char(options.glacier(g));
        %MEAN
         %Intercept
        sweMEAN.(glacier) = repmat(RegressC.(glacier){end,1}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        coeff = RegressC.(glacier).Properties.RowNames(1:end-2);
        for n = 1:length(coeff)
            param               = char(coeff(n));  
            sweT                = topo_full.(glacier).(param)*RegressC.(glacier){n,1};
            sweMEAN.(glacier)   = sweMEAN.(glacier) + sweT;
        end
        
            sweMEAN.(glacier)(sweMEAN.(glacier)<0) = 0;

        
        %MIN
         %Intercept
        sweMIN.(glacier) = repmat(RegressC.(glacier){end,2}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        coeff = RegressC.(glacier).Properties.RowNames(1:end-2);
        for n = 1:length(coeff)
            param               = char(coeff(n));  
            sweT                = topo_full.(glacier).(param)*RegressC.(glacier){n,2};
            sweMIN.(glacier)    = sweMIN.(glacier) + sweT;
        end
            sweMIN.(glacier)(sweMIN.(glacier)<0) = 0;

        %MAX
         %Intercept
        sweMAX.(glacier) = repmat(RegressC.(glacier){end,3}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        coeff = RegressC.(glacier).Properties.RowNames(1:end-2);
        for n = 1:length(coeff)
            param               = char(coeff(n));  
            sweT                = topo_full.(glacier).(param)*RegressC.(glacier){n,3};
            sweMAX.(glacier)    = sweMAX.(glacier) + sweT;
        end
            sweMAX.(glacier)(sweMAX.(glacier)<0) = 0;
        
        %RANGE
        sweRANGE.(glacier)      = sweMAX.(glacier)-sweMIN.(glacier);
end

    clear param sweT *Coeff glacier g n t i coeff 

%% Anova between MLR and BMS

 %Predicted SWE
 clear p
for g = 1:3
    glacier = char(options.glacier(g));
    for r = 2:9
    p(r,g) = anova1([sweBMS(r).(glacier)(:), sweMLR(r).(glacier)(:)],[],'off');   
    end
end

 %Coeffs with density options
for g = 1:3
    glacier = char(options.glacier(g));
    for r = 1:length(boxMLR.(glacier){:,:})
    yBMS = boxBMS.(glacier){r,:}'; 
    yMLR = boxMLR.(glacier){r,:}'; 
    p(r,g) = anova1([yMLR,yBMS],[],'off');   
    end
end    


%% Geotiff histogram
maxH = 30;
minH = -30;

corr = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/diff13-3_corrected.tif');
    corr(corr>maxH) = NaN;    corr(corr<minH) = NaN;

OG = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/diff13-3_original.tif');
    OG(OG>maxH) = NaN;   OG(OG<minH) = NaN;

    clf
hist(OG(:),50); hold on
hist(corr(:),50); 

ho = findobj(gca,'Type','patch');

ho(2,1).FaceColor = rgb('DarkBlue');  ho(2,1).EdgeColor = 'w';
ho(1,1).FaceColor = rgb('MediumTurquoise');  ho(1,1).EdgeColor = 'w';
    ho(1,1).FaceAlpha = 0.6;

    xlabel('Vertical Difference (m)');  ylabel('Frequency');
    legend('Original DEMs','Corrected DEMs', 'Location','NorthWest');
        
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 7.3 7];
    saveFIG('DEMcorrection_hist')
    
display(['original mean difference = ', num2str(nanmean(OG(:)))])
display(['corrected mean difference = ', num2str(nanmean(corr(:)))])

%% Basic linear regression

 %Plot coefficient values
cmap = cbrewer('qual','Dark2',3);
figure(3);  clf
for g = 1:3
    glacier = char(options.glacier(g));
    X = [struct2table(topo_sampled.(glacier)), table(SWE(g).swe,'VariableNames',{'swe'})];
    basicLR.(glacier) = fitlm(X);
     
    coeffs = [basicLR.(glacier).Coefficients{2:8,1}, BMS(9).(glacier){1:7,1}, MLR(9).(glacier){1:7,1}];
    subplot(1,3,g)
    B = bar(coeffs, 'EdgeColor','none');
    for i = 1:3; B(i).FaceColor = cmap(i,:); end
    legend('Basic','BMS','MLR')
    set(gca,'xticklabel',options.topoVars)
end
  
 %Plot RMSE of all fits
figure(2);  clf
for g = 1:3
    glacier = char(options.glacier(g));
    rmseBasic.(glacier) =  sqrt(sum((SWE(g).swe-basicLR.(glacier).Fitted).^2)/numel(SWE(g).swe));
    rmseBMA.(glacier)   =  sqrt(sum(residualsBMS(8).(glacier).^2)/numel(residualsBMS(8).(glacier)));
    rmseMLR.(glacier)   =  sqrt(sum(residualsMLR(8).(glacier).^2)/numel(residualsMLR(8).(glacier)));
    rmseRK.(glacier)    =  sqrt(sum(residualsRK(8).(glacier).^2)/numel(residualsRK(8).(glacier)));
    rmseKRIG.(glacier)  =  sqrt(sum(residualsKRIG(8).(glacier).^2)/numel(residualsKRIG(8).(glacier)));
rmseALL(g,:) = [rmseBasic.(glacier), rmseBMA.(glacier), rmseMLR.(glacier),...
                                        rmseKRIG.(glacier),rmseRK.(glacier)];
end

B = bar(rmseALL, 'EdgeColor','none');
    cmap = cbrewer('qual','Dark2',size(rmseALL,2)+1);
    for i = 1:size(rmseALL,2); B(i).FaceColor = cmap(i+1,:); end
    legend('BasicLR','BMA','MLR', 'SK','RK')
    ylabel('RMSE (m w.e.)'); 
    set(gca,'xticklabel',{'Glacier 4','Glacier 2','Glacier 13'})

