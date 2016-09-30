%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

%% Import Topo Params
[topo, topotext, toporaw] = xlsread('SPOT_TopoParams.xlsx','SPOT_TopoParams','A1:D3936');

div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];
glacier = [4,2,13];

for i = 1:3
name = ['G', num2str(glacier(i))];
    aspect.(name)      = topo(div(i,1):div(i,2),2);
    elevation.(name)   = topo(div(i,1):div(i,2),3);
    slope.(name)       = topo(div(i,1):div(i,2),4);

%Standardizing variables
    aspect.(name)      = (aspect.(name)-mean(aspect.(name)))/std(aspect.(name));
    elevation.(name)   = (elevation.(name)-mean(elevation.(name)))/std(elevation.(name));
    slope.(name)       = (slope.(name)-mean(slope.(name)))/std(slope.(name));
end 

    clear i name topo*

%% Sx import & stepwise MLR

[d300, d300text] = xlsread('sampling_d300_h0.xlsx','sampling_d300_h0','A1:BX3936');
	d300(:,1:3) = [];   d300text(:,1:3) = [];
    d300(:,43) = [];    d300text(:,43) = [];

% [d300text, dirI] = sort(d300text); d300 = d300(:,dirI);

for i = 1:3
    y       = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = d300(div(i,1):div(i,2),:);
   
    [SMLR, ~, ~, inmodel] = stepwisefit(X,y);
    text    = d300text(1,inmodel);
    X1      = X(:,inmodel);
    mlr_Sx  = regress(y,X1); mlr_sort = real(sort(complex(mlr_Sx)));
    
    best            = find(mlr_Sx == mlr_sort(end-1,1));
    winddir.(name)  = text(best);
    Sx.(name)       = X1(:,best);
end
        clear best d300* i inmodel mlr* name SMLR text X* y
%% MLR

for i = 1:3
    y = SWE(i).swe;
    name = ['G', num2str(glacier(i))];
    X = [aspect.(name), elevation.(name), slope.(name), Sx.(name)];

    [mlr.(name), rmse.(name)] = MLRcalval(y, X);
    best = find(rmse.(name)==min(rmse.(name)));
    mlr_best.(name) = mlr.(name)(best,:);   rmse_best.(name) = rmse.(name)(best,1);
end
        clear best i name X y
%% Plots
