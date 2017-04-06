%% Selecting Data from Pattern

load TopoSWE.mat

input.SWE = transectSWE; input.topo_sampled = transectTOPO; 
input.topo_sampled_ns = topo_sampled_ns;

% balanceRK(5).pattern            = 9999;
% balance_residualsRK(5).pattern  = 9999;
% BMAsubset(5).pattern            = 9999;

type = 'centreline';
n = 10;

for S = 1%1:5
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

[ SWEdata, TOPOdata ] = SortNSelect( SWEdata, TOPOdata, n );

%% Plot - locations
    param = 'empty';
    topoParam.G4  = NaN(options.mapsize(1,:));
    topoParam.G2  = NaN(options.mapsize(2,:));
    topoParam.G13 = NaN(options.mapsize(3,:));

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'colour', 'NOmassB')
    saveFIG(['SamplingLocation_', subset,num2str(S)])%,'_numpeople',num2str(T)])

%% Linear regression

subsetLR.(type) =  LinearRegression( SWEdata, TOPOdata, topo_full );
    
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
 
