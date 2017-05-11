%% Full Data
 for d = 2:9
load TopoSWE.mat
clear SWE Density

options.DensitySWE     = d;
run MeasurementLocations.m  %This program determines the easting and northing of transect measurements
run Import_Density.m        %Imports snow density values
run Import_Transect.m       %Imports transect snow depth and measurement location data
run Import_Zigzag.m         %Imports zigzag snow depth and measurement location data
run Import_SWE.m            %Converts to SWE and condences data

den = options.DenOpt{d-1};
for g = 1:3; 
    glacier = options.glacier{g};  fullSWE.(den).(glacier) = SWE(g); 
    fullSWE.(den).input.(glacier) = [fullSWE.(den).(glacier).swe, fullSWE.(den).(glacier).utm(:,1:2), fullSWE.(den).(glacier).cellN];
end
 end
 
for d = 1:8
    den = options.DenOpt{d};
    display(den)
[ tempswe.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled); 

% Linear regression

  fullLR.(den) =  LinearRegression( tempswe.(den), TOPOdata, topo_full );
    
% Simple kriging

  fullSK.(den) =  KrigingR_G( tempswe.(den) );

% Regression Kriging
 
  fullRK.(den) =  RegressionKriging( tempswe.(den), TOPOdata, topo_full, SWE );
end


%% Selecting Data from Pattern

%     subsetLR(length(n)+1).(type)      = 9999;
%     subsetSK(length(n)+1).(type)      = 9999;
%     subsetRK(length(n)+1).(type)      = 9999;
    subset = 'pattern';
    
load TopoSWE.mat
for t = 1%:5
    
if     t == 1; type = 'centreline';          n = 5:5:50;     clt = 1;
elseif t == 2; type = 'CentreTransect4';     n = 10:10:100;  clt = 2;
elseif t == 3; type = 'CentreTransect3';     n = 10:10:100;  clt = 3;
elseif t == 4; type = 'hourglass';           n = 10:10:100;  clt = 4;
elseif t == 5; type = 'hourglassCircle';     n = 10:10:100;  clt = 5;
end


for c = 1:length(n)
    
    for d = 1:8
den = options.DenOpt{d};

input.SWE = fullSWE.(den); input.topo_sampled = topo_sampled; 
input.topo_sampled_ns = topo_sampled_ns;

     display([type, ' n=',num2str(n(c)),' ',den])
%for S = 1%1:5
%     for T = 1:3
 %Pattern
% subset           = 'pattern';
%     S = 2;
% clt       = S;

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


[ subsetSWE(c).(type).(den), TOPOdata ] = DataSubset( subset, clt, input );

[ subsetSWE(c).(type).(den), TOPOdata ] = ObsInCell(subsetSWE(c).(type).(den), TOPOdata); 

[ subsetSWE(c).(type).(den), TOPOdata ] = SortNSelect( subsetSWE(c).(type).(den), TOPOdata, n(c) );

    % Correct the centreline values for invertable matrix when only centreline
    if strcmp(type,'centreline'); TOPOdata.G13.centreD = repmat(0.001, n(c), 1); end
  
    % Add Accumulation area points  
accumulation = 'false';
    if strcmp(accumulation, 'true')
       [sweA, topoA] = DataSubset( subset, 'accum', input );
       for g = 1:3; glacier = options.glacier{g};
           subsetSWE(c).(type).(den).(glacier) = [subsetSWE(c).(type).(den).(glacier); sweA.(glacier)];
           TOPOdata.(glacier).elevation = [TOPOdata.(glacier).elevation;topoA.(glacier).elevation];
           TOPOdata.(glacier).centreD   = [TOPOdata.(glacier).centreD;  topoA.(glacier).centreD];
           TOPOdata.(glacier).aspect    = [TOPOdata.(glacier).aspect;   topoA.(glacier).aspect];
           TOPOdata.(glacier).slope     = [TOPOdata.(glacier).slope;    topoA.(glacier).slope];
           TOPOdata.(glacier).northness = [TOPOdata.(glacier).northness;topoA.(glacier).northness];
           TOPOdata.(glacier).curvature = [TOPOdata.(glacier).curvature;topoA.(glacier).curvature];
           TOPOdata.(glacier).Sx        = [TOPOdata.(glacier).Sx;       topoA.(glacier).Sx];
       end
    end
    
% Plot - locations
%     topoParam.G4  = NaN(options.mapsize(1,:));    topoParam.G2  = NaN(options.mapsize(2,:));    topoParam.G13 = NaN(options.mapsize(3,:));
% PlotTopoParameter(topoParam,'none', 'SWE (m w.e.)', SWEdata(c).(type).(den), 'colour', 'NOmassB')
%     saveFIG(['SamplingLocation_', subset,num2str(S)])%,'_numpeople',num2str(T)])

% Linear regression

  subsetLR(c).(type).(den) =  LinearRegression( subsetSWE(c).(type).(den), TOPOdata, topo_full );
    
% Simple kriging
 
  subsetSK(c).(type).(den) =  KrigingR_G( subsetSWE(c).(type).(den) );

% Regression Kriging
 
  subsetRK(c).(type).(den) =  RegressionKriging( subsetSWE(c).(type).(den), TOPOdata, topo_full, SWE );
    
    end
end
end


%% PLOT -> map estimate

subs = fieldnames(subsetRK);
for s = 1:length(subs)
for c = 1:10;
            %type = 'centreline';          n = 10:5:50;    
            %type = 'CentreTransect4';     n = 10:10:100;  
            %type = 'CentreTransect3';     n = 10:10:100; 
            %type = 'hourglass';           n = 10:10:100; 
            %type = 'hourglassCircle';     n = 10:10:100; 
            type = subs{s};
            den = options.DenOpt{7};

% figure(1); PlotTopoParameter(subsetLR(c).(type).(den),type, 'SWE (m w.e.)', subsetSWE(c).(type).(den), 'black', 'massB')
%      title('Linear Regression')
%      saveFIG(['MapSubset_LR',type,'_n',num2str(n(c)),den])

figure(2); PlotTopoParameter(subsetSK(c).(type).(den),type, 'SWE (m w.e.)', subsetSWE(c).(type).(den), 'black', 'massB')
     title('Simple Kriging')
     saveFIG(['MapSubset_SK',type,'_n',num2str(n(c)),den])

% figure(3); PlotTopoParameter(subsetRK(c).(type).(den),type, 'SWE (m w.e.)', subsetSWE(c).(type).(den), 'black', 'massB')
%      title('Regression Kriging')
%      saveFIG(['MapSubset_RK',type,'_n',num2str(n(c)),den])
end
end
%% PLOT -> sample size and density variation

for g = 1:3
    glacier = options.glacier{g};
    for d = 1:8
        den = options.DenOpt{d};
        meanLR.(glacier)(d) = nanmean(fullLR.(den).(glacier)(:));
        meanSK.(glacier)(d) = nanmean(fullSK.(den).(glacier)(:));
        meanRK.(glacier)(d) = nanmean(fullRK.(den).(glacier)(:));
    end
end

for g = 1:3
    glacier = char(options.glacier(g));  
    for i = 2:9
    stackSWE.MLR.(glacier)(:,:,i-1)   = sweMLR(i).(glacier);
    stackSWE.BMS.(glacier)(:,:,i-1)   = sweBMS(i).(glacier);
    stackSWE.SK.(glacier)(:,:,i-1)    = sweKRIG(i).(glacier).pred;
    stackSWE.RK.(glacier)(:,:,i-1)    = sweRK(i).(glacier);
    end
end    

%type = 'centreline';    n = 10:5:50;
%load Topo_Regress_Krig.mat
%load Subset.mat

subs = fieldnames(subsetRK);
for s = 1:length(subs)
    type = subs{s};
%Linear Regression
figure(5); clf
for g = 1:3
    glacier = options.glacier{g};
    meanLR.(glacier) = nanmean(stackSWE.MLR.(glacier)(:));
    for c = 1:length(n)
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = nanmean(subsetLR(c).(type).(den).(glacier)(:));
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[meanLR.(glacier), meanLR.(glacier)],'k--')
        xlabel('Sample size'); ylabel('Mean SWE')
        title(['LR ',type, ' ',glacier])
        columnlegend(3,[DenOpt, {'All'}],'Location','NorthEast');
        ylim([0 0.8])
end
    saveFIG(['SubsetMeanSWE_samplesizeNdensity_LR',type],12)

 %Simple Kriging
figure(5); clf
for g = 1:3
    glacier = options.glacier{g};
    meanSK.(glacier) = nanmean(stackSWE.SK.(glacier)(:));
    for c = 1:length(n)
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = nanmean(subsetSK(c).(type).(den).(glacier)(:));
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([0, max(n)],[meanSK.(glacier), meanSK.(glacier)],'k--')
        xlabel('Sample size'); ylabel('Mean SWE')
        title(['SK ',type, ' ',glacier])
        columnlegend(3,[DenOpt, {'All'}],'Location','NorthEast');
        ylim([0 0.8]); xlim([0, max(n)])
end
    saveFIG(['SubsetMeanSWE_samplesizeNdensity_SK',type],12)

 %Regression Kriging
figure(5); clf
for g = 1:3
    glacier = options.glacier{g};
    meanRK.(glacier) = nanmean(stackSWE.RK.(glacier)(:));
    for c = 1:length(n)
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = nanmean(subsetRK(c).(type).(den).(glacier)(:));
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[meanRK.(glacier), meanRK.(glacier)],'k--')
        xlabel('Sample size'); ylabel('Mean SWE')
        title(['RK ',type, ' ',glacier])
        columnlegend(3,[DenOpt, {'All'}],'Location','NorthEast');
        ylim([0 0.8])
end
    saveFIG(['SubsetMeanSWE_samplesizeNdensity_RK',type],12)
end

%% PLOT -> elevation coefficient G2 and G13
figure(7); clf
type = 'centreline';    n = 10:5:50;
for g = 2:3
    glacier = options.glacier{g};
for c = 1:length(n)
for d = 1:8 
   den = DenOpt{d};
   stackC.(glacier)(d,c) = (subsetLR(c).(type).(den).coeff{1,g});
end
end

    subplot(1,2,g-1)
    plot(n,stackC.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[mean(boxMLR.(glacier){1,:}), mean(boxMLR.(glacier){1,:})],'k--')
        xlabel('Sample size'); ylabel('Elevation Regression Coefficient')
        title(['LR Centreline ',glacier])
        legend([DenOpt, {'All'}])
end

%% PLOT -> all coefficients for one sample size
        
            c = 10;
            type = 'Ahourglass';    n = 10:10:100;
clear stackC
for g = 1:3
    glacier = options.glacier{g};
for d = 1:8 
   den = DenOpt{d};
   stackC.(glacier)(d,:) = (subsetLR(c).(type).(den).coeff{1:7,g});
end

figure(8);    subplot(1,3,g)
    boxplot(stackC.(glacier), 'Colors', options.RGB(g,:)); hold on
        plot([0,8],[0 0],'k--')
        ylabel('Regression Coefficient')
        title(['LR ',type,' ',glacier])
        set(gca,'XTickLabels',options.topoVars)
end            
     %saveFIG(['MapSubset_RK',type,'_n',num2str(n(c)),den])

     
%% Residuals of subsets

    subset = 'pattern';
% load TopoSWE.mat
% load Subset.mat    

subs = fieldnames(subsetRK);
for s = 1:length(subs)
type = subs{s};
for c = 1:length(n)
    
    for d = 1:8
    den = options.DenOpt{d};

    SWEpredLR(c).(type).(den) = SampledCell( subsetLR(c).(type).(den) );
    SWEpredSK(c).(type).(den) = SampledCell( subsetSK(c).(type).(den) );
    SWEpredRK(c).(type).(den) = SampledCell( subsetRK(c).(type).(den) );
  
    for g = 1:3
    glacier = options.glacier{g};
    subsetResLR(c).(type).(den).(glacier) = SWE(g).swe-SWEpredLR(c).(type).(den).(glacier);
    subsetResSK(c).(type).(den).(glacier) = SWE(g).swe-SWEpredSK(c).(type).(den).(glacier);
    subsetResRK(c).(type).(den).(glacier) = SWE(g).swe-SWEpredRK(c).(type).(den).(glacier);
    
    subsetRmseLR(c).(type).(den).(glacier) = sqrt(mean((SWE(g).swe-SWEpredLR(c).(type).(den).(glacier)).^2));
    subsetRmseSK(c).(type).(den).(glacier) = sqrt(mean((SWE(g).swe-SWEpredSK(c).(type).(den).(glacier)).^2));
    subsetRmseRK(c).(type).(den).(glacier) = sqrt(mean((SWE(g).swe-SWEpredRK(c).(type).(den).(glacier)).^2));

    end

    end
end
end

%% PLOT -> RMSE
load TopoBMS_MLR.mat SWE
run OPTIONS
 %RMSE of full data
    for d = 1:8
    den = options.DenOpt{d};

    resLR(d) = SampledCell( fullLR.(den) );
    resSK(d) = SampledCell( fullSK.(den) );
    resRK(d) = SampledCell( fullRK.(den) );
  
    for g = 1:3
    glacier = options.glacier{g};
    
    rmseLR(g,d) = sqrt(mean((SWE(g).swe-resLR(d).(glacier)).^2));
    rmseSK(g,d) = sqrt(mean((SWE(g).swe-resSK(d).(glacier)).^2));
    rmseRK(g,d) = sqrt(mean((SWE(g).swe-resRK(d).(glacier)).^2));
    end
    end
    
subs = fieldnames(subsetRK);
for s = 1:length(subs)
    type = subs{s};
%Linear Regression
figure(4); clf
for g = 1:3
    glacier = options.glacier{g};
    for c = 1:length(n)
for d = 1:8 
   den = options.DenOpt{d};
   stack.(glacier)(d,c) = subsetRmseLR(c).(type).(den).(glacier);
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[mean(rmseLR(g,:)), mean(rmseLR(g,:))],'k--')
        xlabel('Sample size'); ylabel('RMSE (m w.e.)')
        title(['LR ',type, ' ',glacier])
        legend([DenOpt, {'All'}],'Location','best')
        %ylim([0.05 0.2])
        set(gca,'YTick',(0:0.01:1))
end
    saveFIG(['SubsetRMSE_samplesizeNdensity_LR',type],18,'3G')

 %Simple Kriging
figure(5); clf
for g = 1:3
    glacier = options.glacier{g};
    for c = 1:length(n)
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = subsetRmseSK(c).(type).(den).(glacier);
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[mean(rmseSK(g,:)), mean(rmseSK(g,:))],'k--')
        xlabel('Sample size'); ylabel('RMSE (m w.e.)')
        title(['SK ',type, ' ',glacier])
        legend([DenOpt, {'All'}],'Location','best')
        %ylim([0.05 0.2])
        set(gca,'YTick',(0:0.01:1))
end
    saveFIG(['SubsetRMSE_samplesizeNdensity_SK',type],18,'3G')

 %Regression Kriging
figure(6); clf
for g = 1:3
    glacier = options.glacier{g};
    for c = 1:length(n)
for d = 1:8 
   den = DenOpt{d};
   stack.(glacier)(d,c) = subsetRmseRK(c).(type).(den).(glacier);
end
    end
    subplot(1,3,g)
    plot(n,stack.(glacier),'LineWidth',2); hold on
    plot([min(n), max(n)],[mean(rmseRK(g,:)), mean(rmseRK(g,:))],'k--')
        xlabel('Sample size'); ylabel('RMSE (m w.e.)')
        title(['RK ',type, ' ',glacier])
        legend([DenOpt, {'All'}],'Location','best')
        %ylim([0.05 0.2])
        set(gca,'YTick',(0:0.01:1))
end
    saveFIG(['SubsetRMSE_samplesizeNdensity_RK',type],18,'3G')
end

%% PLOT -> compare sampling designs and interpolation methods over sample size
for g = 1:3
    glacier = options.glacier{g};
    for d = 1:8
        den = options.DenOpt{d};
        meanLR(g,d) = nanmean(fullLR.(den).(glacier)(:));
        meanSK(g,d) = nanmean(fullSK.(den).(glacier)(:));
        meanRK(g,d) = nanmean(fullRK.(den).(glacier)(:));
    end
end
%RMSE of full data
    %SWE = ObsInCell(SWE, topo_sampled);
for d = 1:8
den = options.DenOpt{d};

resLR(d) = SampledCell( fullLR.(den) );
resSK(d) = SampledCell( fullSK.(den) );
resRK(d) = SampledCell( fullRK.(den) );
    for g = 1:3
    glacier = options.glacier{g};

    rmseLR(g,d) = sqrt(mean((SWE(g).swe-resLR(d).(glacier)).^2));
    rmseSK(g,d) = sqrt(mean((SWE(g).swe-resSK(d).(glacier)).^2));
    rmseRK(g,d) = sqrt(mean((SWE(g).swe-resRK(d).(glacier)).^2));
    end
end

figure(1); clf
    subsets = fieldnames(subsetRK); subsets = sort(subsets);
    subsets = subsets([3,1,2,4,5,8,6,7,9,10]);
    pp = 1;    
for r = 1:3
    Xlab = 'RMSE (m w.e.)';
    if     r ==1; D = subsetRmseLR; interp = 'LR'; rmseKK = rmseLR;
    elseif r ==2; D = subsetRmseSK; interp = 'SK'; rmseKK = rmseSK;
    elseif r ==3; D = subsetRmseRK; interp = 'RK'; rmseKK = rmseRK;
    end
%     Xlab = 'SWE (m w.e.)';
%     if     r ==1; D = subsetLR; interp = 'LR'; rmseKK = meanLR;
%     elseif r ==2; D = subsetSK; interp = 'SK'; rmseKK = meanSK;
%     elseif r ==3; D = subsetRK; interp = 'RK'; rmseKK = meanRK;
%     end
for g = 1:3
     glacier = options.glacier{g};
 for s = 1:length(subsets)
     sub = subsets{s};
    for n = 1:9
        for d = 1:8
        den = options.DenOpt{d};
        T(d)= nanmean(D(n).(sub).(den).(glacier)(:));
        end
        rmseS.(glacier)(n,s)= mean(T);
    end
 end
 
     cols = [44, 117, 234;   14, 140, 9;    140, 9, 140;    18, 200, 204;    140, 74, 9];
        cols = cols/255;
     n = 10:10:90;
 subplot(3,3,pp)
    plot(10:5:50,rmseS.(glacier)(:,1),'Color',cols(1,:)); hold on
 for l = 2:5
    plot(n,rmseS.(glacier)(:,l),'Color',cols(l,:)); hold on
 end
     plot(10:5:50,rmseS.(glacier)(:,6),'--','Color',cols(1,:)); hold on
for l = size(rmseS.(glacier),2)/2+2:size(rmseS.(glacier),2)
    plot(n,rmseS.(glacier)(:,l),'--','Color',cols(l-size(rmseS.(glacier),2)/2,:))
end
        ylabel([interp, ' ',Xlab]); xlabel('Sample Size')
        if pp<=3; title(glacier); end
    plot([0, max(n)],[mean(rmseKK(g,:)), mean(rmseKK(g,:))],'k--')    
if strcmp(Xlab,'RMSE (m w.e.)')
        if      g==1; ylim([0.11 0.23])
        elseif  g==2; ylim([0.065 0.14]);
        elseif  g==3; ylim([0.065 0.14]);
        end
        set(gca,'YTick',(0:0.01:1))
elseif strcmp(Xlab,'SWE (m w.e.)')
         ylim([0.2 0.7]);
end
xlim([0 max(n)])
pp = pp+1;
end
end
    legend(subsets);
    saveFIG(['SubsetInterpSizeCompile_',Xlab(1:3)],14)
    
%% PLOT -> sampling designs
    subset = 'pattern';
for t = 1:5
    
if     t == 1; type = 'centreline';          n = 10:5:50;    clt = 1;
elseif t == 2; type = 'CentreTransect4';     n = 10:10:100;  clt = 2;
elseif t == 3; type = 'CentreTransect3';     n = 10:10:100;  clt = 3;
elseif t == 4; type = 'hourglass';           n = 10:10:100;  clt = 4;
elseif t == 5; type = 'hourglassCircle';     n = 10:10:100;  clt = 5;
end


load TopoSWE.mat
input.SWE = fullSWE.(den); input.topo_sampled = topo_sampled; 
input.topo_sampled_ns = topo_sampled_ns;


[ Design.(type), TOPOdata ] = DataSubset( subset, clt, input );

 
    % Add Accumulation area points  
accumulation = 'true';
    if strcmp(accumulation, 'true')
       [sweA, topoA] = DataSubset( subset, 'accum', input );
       for g = 1:3; glacier = options.glacier{g};
           Design.(type).(glacier) = [ Design.(type).(glacier); sweA.(glacier)];       end
    end
    
% Plot - locations
     topoParam.G4  = NaN(options.mapsize(1,:));    topoParam.G2  = NaN(options.mapsize(2,:));    topoParam.G13 = NaN(options.mapsize(3,:));
 PlotTopoParameter(topoParam,'none', 'SWE (m w.e.)',  Design.(type), 'colour', 'NOmassB')
     saveFIG(['SamplingDesign_A', type])
    
end

%% Random subsets

input.SWE               = SWE;  
input.topo_sampled      = topo_sampled; 
input.topo_sampled_ns   = topo_sampled_ns;
n                       = 10:5:60;

for i = 1:length(n)
    for x = 1:30
        [ subsetSWE, TOPOdata ] = DataSubset( 'random', n(i), input );

        [ testRK ] = RegressionKriging( subsetSWE, TOPOdata, topo_full, SWE );
        for g = 1:3;
        glacier = char(options.glacier(g));
        Ttemp = KrigingR(subsetSWE.(glacier)(:,1), subsetSWE.(glacier)(:,2:3), glacier);
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
 
