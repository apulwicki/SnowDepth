%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%


%% Import Topo Params
%Run MAIN
MAIN

%Import topos
[topo, ~, ~] = xlsread('SPOT_TopoParams.xlsx','sampling_TopoParams','A1:G3936');

div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];
glacier = [4,2,13];

for i = 1:3
name = ['G', num2str(glacier(i))];
    aspect           = topo(div(i,1):div(i,2),2);
    northness        = topo(div(i,1):div(i,2),3);
    profileCurve     = topo(div(i,1):div(i,2),4);
    tangentCurve     = topo(div(i,1):div(i,2),5);
    slope            = topo(div(i,1):div(i,2),6);
    elevation        = topo(div(i,1):div(i,2),7);
    
    topo_sampled.(name) = struct('aspect',aspect, 'elevation',elevation,...
        'northness',northness, 'profileCurve',profileCurve, 'slope',slope,...
        'tangentCurve',tangentCurve);

end

% Sx import & stepwise MLR
[d300, d300text] = xlsread('d300_h0.xlsx','sampling_d300_h0','B1:BU3936');
[d200, d200text] = xlsread('d200_h0.xlsx','sampling_d200_h0','A1:BT3936');
[d100, d100text] = xlsread('d100_h0.xlsx','sampling_d100_h0','A1:BT3936');
    d300text = strcat('d300',d300text);
    d200text = strcat('d200',d200text);
    d100text = strcat('d100',d100text);

for i = 1:3
    y       = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = [d100(div(i,1):div(i,2),:), d200(div(i,1):div(i,2),:), d300(div(i,1):div(i,2),:)];
    text    = [d100text, d200text, d300text];
   
    [~, ~, ~, inmodel] = stepwisefit(X,y);
    text    = text(1,inmodel);
    X1      = [ones(length(X),1), X(:,inmodel)];
    mlr_Sx     = regress(y,X1); mlr_sort = real(sort(complex(mlr_Sx)));
    
    best            = find(mlr_Sx == mlr_sort(end-1,1));
    winddir.(name)  = text(best);
    Sx              = X1(:,best);
    topo_sampled.(name).Sx = Sx;
end
        clear best d300* d200* d100* i inmodel mlr* name text X* y

%Distance from centreline import
run DistanceFromCentreline.m

for i = 1:3
   name     = char(glacier(i));
   topo_sampled.(name).centreD = min_dist.(name);
end
            
%Standardizing variables
params = fieldnames(topo_sampled.G4);
for i = 1:3
name= char(glacier(i));
    for t = 1:length(params)
    field = char(params(t));
    
    topo_sampled.(name).(field) = (topo_sampled.(name).(field)-...
        mean(topo_sampled.(name).(field)))/std(topo_sampled.(name).(field));
    end
end

    clear i name topo field j params div glacier
    clear aspect* elev* north* profil* slope* Sx* tangent*

%% MLR - Topo Regression

glacier = {'G4','G2','G13'};
%remove aspect
for i = 1:3
   name     = char(glacier(i));
   topo_sampled.(name) = rmfield(topo_sampled.(name),'aspect');
end

for t = 2:9
options.DensitySWE = t;
run MeasurementLocations.m  %This program determines the easting and northing of transect measurements
run Import_Density.m %Imports snow density values
run Import_Transect.m %Imports transect snow depth and measurement location data
run Import_Zigzag.m %Imports zigzag snow depth and measurement location data
run Import_SWE.m %Converts to SWE and condences data

    for i = 1:3
        glacier = [4,2,13];
        y       = SWE(i).swe;
        name    = ['G', num2str(glacier(i))]; 
            display(['option = ',num2str(t), ', glacier = ',name]);
        X       = topo_sampled.(name);

        %No zigzag
        NoZZ    = SWE(i).pattern~='ZZ';
        y = y(NoZZ);    
        f = fieldnames(X);
        for d = 1:length(f)
            param = char(f(d));     X.(param) = X.(param)(NoZZ);
        end
        
        [mlr(t).(name), rmse(t).(name), lm(t).(name)] = MLRcalval(y, X);
        mlr(t).(name).Properties.VariableNames = strcat(mlr(t).(name).Properties.VariableNames, num2str(t));
    end
end
        clear best i name X y t glacier

%% Export all values
G4_mlrDensity = []; G2_mlrDensity = []; G13_mlrDensity = [];
for i = 2:9
G4_mlrDensity  = [G4_mlrDensity, mlr(i).G4(:,3)];
G2_mlrDensity  = [G2_mlrDensity, mlr(i).G2(:,3)];
G13_mlrDensity = [G13_mlrDensity, mlr(i).G13(:,3)];
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
%% Plots - MLR

% Actual vs fitted data
figure
    axis([0 1.2 0 1.2]);
line = refline(1,0);
    line.Color = 'k'; line.LineStyle = '--'; hold on
RGB = [9, 132, 103; 224, 187, 2; 130, 75, 135]/255;

option = 8;

for i = 1:3
    y  = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = topo_sampled.(name);
    %No zigzag
    NoZZ    = SWE(i).pattern~='ZZ';
    y = y(NoZZ);    
    ff = fieldnames(X);
    for d = 1:length(ff)
        param = char(ff(d));     X.(param) = X.(param)(NoZZ);
    end
        params = fieldnames(X); 
        M = X.(char(params(1)));
        for j = 2:length(params)
            field = char(params(j));
            M = [M, X.(field)];
        end
fitted_swe = sum(repmat(table2array(mlr(option).(name)(:,1))',length(M),1).*...
                        [ones(length(M),1),M],2);

    plot(y, fitted_swe, '.', 'Color', RGB(i,:),'MarkerSize',13); hold on

    [f.(name), g.(name)] = fit(y, fitted_swe,'poly1');
    p = plot(f.(name)); hold on
    set(p,'Color',RGB(i,:)); set(p, 'LineWidth',1.5);    
end
    xlabel('Original SWE (m)'); ylabel('MLR SWE (m)');
    legend('Reference Line','Glacier 4',['R^2=',num2str(round(g.G4.rsquare,2))],...
        'Glacier 2',['R^2=',num2str(round(g.G2.rsquare,2))],...
        'Glacier 13',['R^2=',num2str(round(g.G13.rsquare,2))],'Location','northwest')    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18) 
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 10 9];

filename = 'MLRfit';
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

%% Plots - Box and whisker for coeffs from density options

%Rearrange to compare density options
params = mlr(2).G4.Properties.RowNames;
box.G4 = []; box.G2 = []; box.G13 = [];
for i = 2:9
    box.G4 = [box.G4; table2array(mlr(i).G4(:,3))'];
    box.G2 = [box.G2; table2array(mlr(i).G2(:,3))'];
    box.G13 = [box.G13; table2array(mlr(i).G13(:,3))'];
end
%No intercept option
box.G4 = box.G4(:,2:end); box.G2 = box.G2(:,2:end); box.G13 = box.G13(:,2:end); params = params(2:end);

h = cat(1, reshape(box.G4,[1 size(box.G4)]), reshape(box.G2,[1 size(box.G2)]),...
            reshape(box.G13,[1 size(box.G13)]));

    RGB = [9, 132, 103; 224, 187, 2; 130, 75, 135]/255;
aboxplot(h,'labels',params, ...
    'Colormap',                 RGB,...
    'OutlierMarkerSize',        10,...
    'OutlierMarkerEdgeColor',   [0 0 0]); % Advanced box plot
        legend('Glacier 4','Glacier 2','Glacier 13'); % Add a legend
        ylabel('% Variance Explained')
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 10 10];

filename = 'Coeffs_DensityOpts';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')

%% AIC values

glacier = {'G4','G2','G13'};

for g = 1:length(glacier)
    G = char(glacier(g));
    for t = 2:9
       aic.(G)(t-1,1)   = t; %density option
       aic.(G)(t-1,2)   = lm(t).(G).ModelCriterion.AIC; %AIC value
    end

    aic.(G)(:,3)        = exp((min(aic.(G)(:,2))-aic.(G)(:,2))/2); %
end
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


%% Regression tree

header = {'aspect','northness', 'profileCurve', ...
                'tangentCurve', 'slope', 'elevation', 'Sx'};
figure;
line = refline(1,0);
    line.Color = 'k'; line.LineStyle = '--'; hold on

for i = 1:3
    y       = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = table(aspect.(name), northness.(name), profileCurve.(name), ...
                tangentCurve.(name), slope.(name), elevation.(name), Sx.(name),...
                'VariableNames',header);

    cal_ind = randperm(length(y),floor(length(y)*3/4));
    val_ind = setdiff(1:length(y),cal_ind);
    
    tree            = fitrtree(X(cal_ind,:),y(cal_ind,1));
    [~,~,~,bestlevel] = cvLoss(tree,'SubTrees','All','TreeSize','min') %minimizing cross-validated loss
    tree = prune(tree,'Level',6);
    %view(tree,'Mode','Graph')
    tree_predict    = predict(tree,X(val_ind,:));
    
    
    plot(y(val_ind),tree_predict,'.'); hold on
end
    
    
    forest = TreeBagger(50, X, y,'OOBPrediction','On','Method','regression');
    view(forest.Trees{1},'Mode','graph')
    oobErrorBaggedEnsemble = oobError(forest);
    plot(oobErrorBaggedEnsemble)
    
    
%% Range of params sampled
header = {'Sx', 'aspect', 'elevation', 'northness', 'profileCurve', 'slope', 'tangentCurve'};
%topo_sampled = [Sx, aspect, elevation, northness, profileCurve, slope, tangentCurve];

%Importing DEMs
files = dir('/home/glaciology1/Documents/Data/GlacierTopos/*.asc');
v = cell(0,0);
for i = 1:length(files)
    A = importdata(['/home/glaciology1/Documents/Data/GlacierTopos/', files(i).name]);
    A.data(A.data==-9999) = NaN;   A.data(1:6) = [];
    v(i,1) = cellstr(files(i).name(1:end-4)); eval([v{i,1} '= A.data;']);
end
    elevation_G04(elevation_G04==0) = NaN; elevation_G02(elevation_G02==0) = NaN; elevation_G13(elevation_G13==0) = NaN;
    northness_G04(northness_G04<-3e+38) = NaN; northness_G02(northness_G02<-3e+38) = NaN; northness_G13(northness_G13<-3e+38) = NaN;
    
for i = 1:length(v)/3
    param = char(v(i*3)); param = param(1:end-4);
    eval(['topo_full.G4.(param) = ',v{3*i-1,1},';']);
    eval(['topo_full.G2.(param) = ',v{3*i-2,1}]);
    eval(['topo_full.G13.(param) = ',v{3*i,1}]);
end
    clear aspect* elev* north* profil* slope* Sx* tangent* files A i param v
   
glacier = {'G4','G2','G13'};
for r = 1:7
figure
    for i = 1:3
        name    = char(glacier(i)); param = char(header(r));
        [N, edges] = histcounts(topo_sampled.(name).(param),10); 
        [Nall, edgesall] = histcounts(topo_full.(name).(param),10); 
        subplot(1,3,i)
            plot((edges(1:end-1)+edges(2:end))/2,N); hold on %sampled values
            plot((edgesall(1:end-1)+edgesall(2:end))/2,Nall)
            xlabel([header(r), ' ', name]); ylabel('Freq.')
            legend('Sampled','Full range')
    end
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 12 6];
filename = ['SampledRangeTopo_',header{r}];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
end 

    clear v i r name header glacier


