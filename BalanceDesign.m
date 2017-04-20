%% Selecting Data from Pattern
    
type = 'hourglass';
n = 10:10:100;
DenOpt = {'S1','F1','S2','F2','S3','F3','S4','F4'};
    subsetLR(length(n)+1).(type)      = 9999;
    subsetSK(length(n)+1).(type)      = 9999;
    subsetRK(length(n)+1).(type)      = 9999;


for c = 1:length(n)
    
    for d = 2:9
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

%for S = 1%1:5
%     for T = 1:3
 %Pattern
subset           = 'pattern';
    S = 6;
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

[ SWEdata(c).(type).(den), TOPOdata ] = DataSubset( subset, option, input );

[ SWEdata(c).(type).(den), TOPOdata ] = ObsInCell(SWEdata(c).(type).(den), TOPOdata);

[ SWEdata(c).(type).(den), TOPOdata ] = SortNSelect( SWEdata(c).(type).(den), TOPOdata, n(c) );

    if strcmp(type,'centreline'); TOPOdata.G13.centreD = repmat(0.001, n(c), 1); end
        
% Plot - locations
%     param = 'empty';
%     topoParam.G4  = NaN(options.mapsize(1,:));
%     topoParam.G2  = NaN(options.mapsize(2,:));
%     topoParam.G13 = NaN(options.mapsize(3,:));
% 
% PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata(c).(type).(den), 'colour', 'NOmassB')
%     saveFIG(['SamplingLocation_', subset,num2str(S)])%,'_numpeople',num2str(T)])

% Linear regression

subsetLR(c).(type).(den) =  LinearRegression( SWEdata(c).(type).(den), TOPOdata, topo_full );
    
% Simple kriging

subsetSK(c).(type).(den) =  KrigingR_G( SWEdata(c).(type).(den) );

% Regression Kriging
 
subsetRK(c).(type).(den) =  RegressionKriging( SWEdata(c).(type).(den), TOPOdata, topo_full, SWE );

%     end
    end
end

%% PLOT -> map estimate

        n = 10:5:50;
        DenOpt = {'S1','F1','S2','F2','S3','F3','S4','F4'};
        
            c = 1;
            type = 'centreline';
            den = DenOpt{7};

figure(1); PlotTopoParameter(subsetLR(c).(type).(den),type, 'SWE (m w.e.)', SWEdata(c).(type).(den), 'black', 'massB')
     title('Linear Regression')
     saveFIG(['MapSubset_LR',type,'_n',num2str(n(c)),den])

figure(2); PlotTopoParameter(subsetSK(c).(type).(den),type, 'SWE (m w.e.)', SWEdata(c).(type).(den), 'black', 'massB')
     title('Simple Kriging')
     saveFIG(['MapSubset_SK',type,'_n',num2str(n(c)),den])

figure(3); PlotTopoParameter(subsetRK(c).(type).(den),type, 'SWE (m w.e.)', SWEdata(c).(type).(den), 'black', 'massB')
     title('Regression Kriging')
     saveFIG(['MapSubset_RK',type,'_n',num2str(n(c)),den])

%% PLOT -> sample size and density variation
n = 10:5:50;

load Topo_Regress_Krig.mat
 %Linear Regression
figure(4); clf
for g = 1:3
    glacier = options.glacier{g};
    meanLR.(glacier) = nanmean(stackSWE.MLR.(glacier)(:));
    for c = 1:length(n);
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = nanmean(subsetLR(c).centreline.(den).(glacier)(:));
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[meanLR.(glacier), meanLR.(glacier)],'k--')
        xlabel('Sample size'); ylabel('Mean SWE')
        title(['LR Centreline ',glacier])
        legend([DenOpt, {'All'}])
        %ylim([0 2.45])
end

 %Simple Kriging
figure(5); clf
for g = 1:3
    glacier = options.glacier{g};
    meanSK.(glacier) = nanmean(stackSWE.SK.(glacier)(:));
    for c = 1:length(n);
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = nanmean(subsetSK(c).centreline.(den).(glacier)(:));
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[meanSK.(glacier), meanSK.(glacier)],'k--')
        xlabel('Sample size'); ylabel('Mean SWE')
        title(['SK Centreline ',glacier])
        legend([DenOpt, {'All'}])
        %ylim([0 2.45])
end

 %Regression Kriging
figure(6); clf
for g = 1:3
    glacier = options.glacier{g};
    meanRK.(glacier) = nanmean(stackSWE.RK.(glacier)(:));
    for c = 1:length(n);
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = nanmean(subsetRK(c).centreline.(den).(glacier)(:));
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[meanRK.(glacier), meanRK.(glacier)],'k--')
        xlabel('Sample size'); ylabel('Mean SWE')
        title(['RK Centreline ',glacier])
        legend([DenOpt, {'All'}])
        %ylim([0 2.45])
end

%% PLOT -> elevation coefficient G2 and G13
figure(7); clf
n = 10:5:50;
for g = 2:3
    glacier = options.glacier{g};
for c = 1:length(n);
for d = 1:8 
   den = DenOpt{d};
   stackC.(glacier)(d,c) = (subsetLR(c).centreline.(den).coeff{1,g});
end
end

    subplot(1,2,g-1)
    plot(n,stackC.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[mean(boxMLR.(glacier){1,:}), mean(boxMLR.(glacier){1,:})],'k--')
        xlabel('Sample size'); ylabel('Elevation Regression Coefficient')
        title(['LR Centreline ',glacier])
        legend([DenOpt, {'All'}])
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
 
