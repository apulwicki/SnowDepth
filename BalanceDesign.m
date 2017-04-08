%% Selecting Data from Pattern

type = 'centreline';
n = 10:5:50;
DenOpt = {'S1','F1','S2','F2','S3','F3','S4','F4'};

for c = 1%:length(n)
    
    for d = 2%:9;
load TopoSWE.mat
clear SWE Density

options.DensitySWE     = d;
run MeasurementLocations.m  %This program determines the easting and northing of transect measurements
run Import_Density.m        %Imports snow density values
run Import_Transect.m       %Imports transect snow depth and measurement location data
run Import_Zigzag.m         %Imports zigzag snow depth and measurement location data
run Import_SWE.m            %Converts to SWE and condences data
transectSWE = SWE;  transectTOPO = topo_sampled;

input.SWE = transectSWE; input.topo_sampled = transectTOPO; 
input.topo_sampled_ns = topo_sampled_ns;
den = DenOpt{d-1};

% balanceRK(5).pattern            = 9999;
% balance_residualsRK(5).pattern  = 9999;
% BMAsubset(5).pattern            = 9999;




%for S = 1%1:5
%     for T = 1:3
 %Pattern
subset           = 'pattern';
    S = 1;
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

[ SWEdata, TOPOdata ] = SortNSelect( SWEdata, TOPOdata, n(c) );

%% Plot - locations
%     param = 'empty';
%     topoParam.G4  = NaN(options.mapsize(1,:));
%     topoParam.G2  = NaN(options.mapsize(2,:));
%     topoParam.G13 = NaN(options.mapsize(3,:));
% 
% PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'colour', 'NOmassB')
%     saveFIG(['SamplingLocation_', subset,num2str(S)])%,'_numpeople',num2str(T)])

%% Linear regression

subsetLR(c).(type).(den) =  LinearRegression( SWEdata, TOPOdata, topo_full );
    
%% Simple kriging

subsetSK(c).(type).(den) =  KrigingR_G( SWEdata );

%% Regression Kriging
 
subsetRK(c).(type).(den) =  RegressionKriging( SWEdata, TOPOdata, topo_full, SWE );

%     end
    end
end

%% Plot estimate

    param = 'RK';
    topoParam.G4  = balanceRK(S).(subset).G4;
    topoParam.G2  = balanceRK(S).(subset).G2;
    topoParam.G13 = balanceRK(S).(subset).G13;

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'black', 'massB')
    saveFIG(['Map_',param, subset,num2str(S)])%,'_numpeople',num2str(T)])
    
%% Random subsets

input.SWE               = SWE;  
input.topo_sampled      = topo_sampled; 
input.topo_sampled_ns   = topo_sampled_ns;
n                       = 10:5:60;

for i = 1:length(n)
    for x = 1:30
        [ SWEdata, TOPOdata ] = DataSubset( 'random', n(i), input );

        [ testRK ] = RegressionKriging( SWEdata, TOPOdata, topo_full, SWE );
        for g = 1:3;
        glacier = char(options.glacier(g));
        Ttemp = KrigingR(SWEdata.(glacier)(:,1), SWEdata.(glacier)(:,2:3), glacier);
        testKRIG.(glacier) = Ttemp.pred;
        end
        
        locRK       = SampledCell(testRK);
        outRK       = RMSE( locRK, SWE );
        tempRK(x,:) = struct2table(outRK);
        
        locKRIG     = SampledCell(testKRIG);
        outKRIG     = RMSE( locKRIG, SWE );
        tempKRIG(x,:) = struct2table(outKRIG);
    end
rmseRK2(i,:)   = mean(tempRK{:,:});
rmseKRIG2(i,:) = mean(tempKRIG{:,:});
end


 figure(1); clf
for g = 1:3;
glacier = char(options.glacier(g));

subplot(1,3,g) 
    plot(n, rmse(:,g),'Color',options.RGB(g,:),'LineWidth',3); hold on; 
%     rl = refline; rl.Color = 'k'; rl.LineStyle = '--';
end

 figure(1); clf
for g = 1:3;
glacier = char(options.glacier(g));

subplot(1,3,g) 
    plot(n, rmseKRIG(:,g),'LineWidth',3); hold on; 
    plot(n, rmseRK(:,g),'LineWidth',3); hold on; 
    legend('SK','RK')
    xlabel('Sample size'); ylabel('RMSE (m w.e.)')
end
 
