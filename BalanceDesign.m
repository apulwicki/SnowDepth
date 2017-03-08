%% Selecting Data from Pattern

load TopoSWE.mat

input.SWE = SWE; input.topo_sampled = topo_sampled; 
input.topo_sampled_ns = topo_sampled_ns;

balanceRK(5).pattern            = 9999;
balance_residualsRK(5).pattern  = 9999;
BMAsubset(5).pattern            = 9999;

n = 1;

for S = 1:5
%     for T = 1:3
 %Pattern
subset           = 'pattern';
option.clt       = S;

%  %Measurement density
% subset           = 'density';
% option.density   = S;
% option.people    = T;
% 
%  %Topographic parameter
% subset           = 'topoparam';
% option.topo      = 'elevation';
% option.lessgreat = 'less';
% option.value     = 2200;

[ SWEdata, TOPOdata ] = DataSubset( subset, option, input );

[ SWEdata, TOPOdata ] = ObsInCell(SWEdata, TOPOdata);
%% Plot - locations
    param = 'empty';
    topoParam.G4  = NaN(options.mapsize(1,:));
    topoParam.G2  = NaN(options.mapsize(2,:));
    topoParam.G13 = NaN(options.mapsize(3,:));

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'colour', 'NOmassB')
    saveFIG(['SamplingLocation_', subset,num2str(S)])%,'_numpeople',num2str(T)])

%% Regression Kriging
 
[ balanceRK(n).(subset), balance_residualsRK(n).(subset), BMAsubset(n).(subset) ] = ...
    RegressionKriging( SWEdata, TOPOdata, topo_full, SWE );

%% Plot estimate

    param = 'RK';
    topoParam.G4  = balanceRK(S).(subset).G4;
    topoParam.G2  = balanceRK(S).(subset).G2;
    topoParam.G13 = balanceRK(S).(subset).G13;

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'black', 'massB')
    saveFIG(['Map_',param, subset,num2str(S)])%,'_numpeople',num2str(T)])
    
    n = n+1;
%     end
end