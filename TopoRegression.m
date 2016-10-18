%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

%Run MAIN
MAIN

%% Import Topo Params
[topo, topotext, toporaw] = xlsread('SPOT_TopoParams.xlsx','sampling_TopoParams','A1:G3936');

div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];
glacier = [4,2,13];

for i = 1:3
name = ['G', num2str(glacier(i))];
    aspect.(name)           = topo(div(i,1):div(i,2),2);
    northness.(name)        = topo(div(i,1):div(i,2),3);
    profileCurve.(name)     = topo(div(i,1):div(i,2),4);
    tangentCurve.(name)     = topo(div(i,1):div(i,2),5);
    slope.(name)            = topo(div(i,1):div(i,2),6);
    elevation.(name)        = topo(div(i,1):div(i,2),7);

%Standardizing variables
    aspect.(name)           = (aspect.(name)-mean(aspect.(name)))/std(aspect.(name));
    northness.(name)        = (northness.(name)-mean(northness.(name)))/std(northness.(name));
    profileCurve.(name)     = (profileCurve.(name)-mean(profileCurve.(name)))/std(profileCurve.(name));
    tangentCurve.(name)     = (tangentCurve.(name)-mean(tangentCurve.(name)))/std(tangentCurve.(name));
    slope.(name)            = (slope.(name)-mean(slope.(name)))/std(slope.(name));
    elevation.(name)        = (elevation.(name)-mean(elevation.(name)))/std(elevation.(name));
end 

    clear i name topo*

%% Sx import & stepwise MLR

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
    Sx.(name)       = X1(:,best);
end
        clear best d300* d200* d100* i inmodel mlr* name text X* y
%% MLR

for i = 1:3
    y       = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = [aspect.(name), northness.(name), profileCurve.(name), ...
                tangentCurve.(name), slope.(name), elevation.(name), Sx.(name)];

    [mlr.(name), rmse.(name)] = MLRcalval(y, X);
    best = find(rmse.(name)==min(rmse.(name)));
    mlr_best.(name) = mlr.(name)(best,:);   rmse_best.(name) = rmse.(name)(best,1); 
end
        clear best i name X y
        
%ANOVA
% for i = 1:3
%     OG_swe  = SWE(i).swe;
%     name    = ['G', num2str(glacier(i))];
%     X       = [aspect.(name), northness.(name), profileCurve.(name), ...
%                 tangentCurve.(name), slope.(name), elevation.(name), Sx.(name)];
% %    fitted_swe = sum(repmat(mlr_best.(name),length(X),1).*[ones(length(X),1),X],2);
% 
% %     display([name, ' ANOVA'])
% %     anova1(OG_swe, fitted_swe)
%     y       = SWE(i).swe;
%     mdl = fitlm(X,y);
%     tbl = anova(mdl)
% end
%     close all
%% Plots - MLR

% Actual vs fitted data
figure
    axis([0 1.2 0 1.2]); box
line = refline(1,0);
    line.Color = 'k'; line.LineStyle = '--'; hold on
RGB = [0 76 153; 0 153 76; 255 127 0]/255;
for i = 1:3
    OG_swe  = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = [aspect.(name), northness.(name), profileCurve.(name), ...
                tangentCurve.(name), slope.(name), elevation.(name), Sx.(name)];
    fitted_swe = sum(repmat(mlr_best.(name),length(X),1).*[ones(length(X),1),X],2);

    plot(OG_swe, fitted_swe, '.', 'Color', RGB(i,:),'MarkerSize',13); hold on

    [f.(name), g.(name)] = fit(OG_swe, fitted_swe,'poly1');
    p = plot(f.(name)); hold on
    set(p,'Color',RGB(i,:)); set(p, 'LineWidth',1.5);    
end
    xlabel('Original SWE (m)'); ylabel('MLR SWE (m)');
    legend('Reference Line','Glacier 4',['R^2=',num2str(round(g.G4.rsquare,2))],...
        'Glacier 2',['R^2=',num2str(round(g.G2.rsquare,2))],...
        'Glacier 13',['R^2=',num2str(round(g.G13.rsquare,2))],'Location','best')    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18) 
    
filename = 'MLRfit';
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

%% ANOVA between topographic params

for i = 1:3
    OG_swe  = SWE(i).swe;
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
header = {'aspect','northness', 'profileCurve', ...
                'tangentCurve', 'slope', 'elevation', 'Sx'};
for i = 1:3
    name    = ['G', num2str(glacier(i))];
    topo.(name) = [aspect.(name), northness.(name), profileCurve.(name), ...
                tangentCurve.(name), slope.(name), elevation.(name), Sx.(name)];
end

for r = 1:7
figure
    for i = 1:3
        name    = ['G', num2str(glacier(i))];
        [N, edges] = histcounts(topo.(name)(:,r));
        subplot(1,3,i)
            plot((edges(1:end-1)+edges(2:end))/2,N) %sampled values
            xlabel(header(r)); ylabel('Freq.')
    end
end    