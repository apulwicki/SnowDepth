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
for i = 1:3
    OG_swe  = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = [aspect.(name), northness.(name), profileCurve.(name), ...
                tangentCurve.(name), slope.(name), elevation.(name), Sx.(name)];
    fitted_swe = sum(repmat(mlr_best.(name),length(X),1).*[ones(length(X),1),X],2);

    display([name, ' ANOVA'])
    anova1(OG_swe, fitted_swe)
end
    close all
%% Plots

% Actual vs fitted data
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
    legend('Glacier 4',['R^2=',num2str(round(g.G4.rsquare,2))],...
        'Glacier 2',['R^2=',num2str(round(g.G2.rsquare,2))],...
        'Glacier 13',['R^2=',num2str(round(g.G13.rsquare,2))],'Location','best')    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',20) 
    
 
    
    
    
    
    
    