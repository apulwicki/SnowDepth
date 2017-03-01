%% Selecting Data from Pattern

input.SWE = SWE; input.topo_sampled = topo_sampled; 
input.topo_sampled_ns = topo_sampled_ns;

 %Pattern
subset           = 'pattern';
option.clt       = 1;

 %Measurement density
subset           = 'density';
option.people    = 2;
option.density   = 3;

 %Topographic parameter
subset           = 'topoparam';
option.topo      = 'elevation';
option.lessgreat = 'less';
option.value     = 2200;

[ SWEdata, TOPOdata, I ] = DataSubset( subset, option, input );

SWEdata = ObsInCell(SWEdata);
%% Plot - locations
    param = 'empty';
    topoParam.G4  = NaN(options.mapsize(1,:));
    topoParam.G2  = NaN(options.mapsize(2,:));
    topoParam.G13 = NaN(options.mapsize(3,:));

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'colour')
    saveFIG('SamplingLocation_subset')

%% Regression Kriging
 
[ balanceRK, balance_residualsRK, BMAsubset ] = ...
    RegressionKriging( SWEdata, TOPOdata, I, topo_full, SWE );

%% Plot estimate

    param = 'RegressionKriging_subset';
    topoParam.G4  = balanceRK.G4;
    topoParam.G2  = balanceRK.G2;
    topoParam.G13 = balanceRK.G13;

    PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'black')
    
     %Save figure
    %saveFIG(['Map_',param])
 