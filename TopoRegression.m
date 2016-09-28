%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

%% Import Topo Params
[topo, topotext, toporaw] = xlsread('SPOT_TopoParams.xlsx','SPOT_TopoParams','A1:D3981');

div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];

for i = 1:3
name = ['g', num2str(i)];
    aspect.(name)      = topo(div(i,1):div(i,2),1);
    elevation.(name)   = topo(div(i,1):div(i,2),2);
    slope.(name)       = topo(div(i,1):div(i,2),3);

%Standardizing variables
% aspect      = (aspect.(name)-mean(aspect.(name)))/std(aspect.(name));
% elevation   = (elevation.(name)-mean(elevation.(name)))/std(elevation.(name));
% slope       = (slope.(name)-mean(slope.(name)))/std(slope.(name));
end 

    clear div i name topo*
%% MLR

for i = 1%:3
    y = SWE(i).swe;
    name = ['g', num2str(i)];
    X = [aspect.(name), elevation.(name), slope.(name)];

    mlr.(name) = MLRcalval(y, X);
end
%% Plots

