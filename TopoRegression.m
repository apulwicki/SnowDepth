%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

%% Import Topo Params
[topo, topotext, toporaw] = xlsread('SPOT_TopoParams.xlsx','SPOT_TopoParams','A1:D3981');


z = pulldataSWE(SWE, 'all','all','all','all');

aspect      = topo(:,1);
elevation   = topo(:,2);
slope       = topo(:,3);

%Standardizing variables
% aspect      = (topo(:,1)-mean(topo(:,1)))/std(topo(:,1));
% elevation   = (topo(:,2)-mean(topo(:,2)))/std(topo(:,2));
% slope       = (topo(:,3)-mean(topo(:,3)))/std(topo(:,3));

%% MLR

% plot(z.swe); hold on; plot(aspect); hold on; plot(elevation); hold on; plot(slope);

y = z.swe;
X = [ones(length(y),1),aspect,elevation,slope];

% applying stepwise regression on Y and X
[a, ~, ~, ~, STATS]=stepwisefit(X(:,2:4),y);
% constant coefficient
a0 = STATS.intercept;

y_regress = a0 + a(1,1)*X(:,1) + a(2,1)*X(:,2);

[f, gof] = fit(y,y_regress,'poly1')
plot(f,y,y_regress)
    title('Stepwise MLR - aspect, elevation, slope'); xlabel('Original SWE (m w.e.)'); xlabel('Regressed SWE (m w.e.)'); 


%% Plots

