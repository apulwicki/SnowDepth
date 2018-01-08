% **********************************************************
% Kriging
% **********************************************************

clear
close all
load TopoBMS_MLR.mat
%% Dice Kriging -> SWE

sweKRIG(9).G4.pred = 9999; 
for g = 2%1:3
        glacier = options.glacier{g};
    for r = 7%2:9
    sweKRIG(r).(glacier) = KrigingR(sweOPT(r).(glacier)(:,1), SWE(g).utm, glacier);
    end
end    
    clear g glacier r

%% 3D plot    
figure(4); clf
for g = 1:3
    r=2;
    glacier = char(options.glacier(g));
subplot(1,3,g)
    surf(sweKRIG(2).(glacier).pred, 'FaceAlpha',0.7,'LineStyle','none'); colorbar; hold on
%     surf(sweKRIG(2).(glacier).upper95, 'FaceColor','r','FaceAlpha',0.3,'LineStyle','none'); colorbar; hold on
%     surf(sweKRIG(2).(glacier).lower95, 'FaceColor','r','FaceAlpha',0.3,'LineStyle','none'); colorbar; hold on
    plot3(options.E.(glacier),options.N.(glacier),sweOPT(r).(glacier)(:,1),'k.')
end
    
%% Plotting -> SWE
    param = 'pred';
    r = 7;
    topoParam.G4  = sweKRIG(r).G4.(param);
    topoParam.G2  = sweKRIG(r).G2.(param);
    topoParam.G13 = sweKRIG(r).G13.(param);
figure(3)
    PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWE, 'black','massB')

    %saveFIG('sweKriged','3G')
    
%% REGRESSION KRIGING

test = KrigingR(residualsBMS(4).G4, SWE(1).utm, 'G4');

% Dice Kriging -> residuals
res_bmaKRIG(9).G4.pred = 9999;
for g = 1:3
        glacier = char(options.glacier(g));
    for r = 4%2:9
    res_bmaKRIG(r).(glacier) = KrigingR(residualsBMS(r).(glacier), SWE(g).utm, glacier);
    end
end    
    clear g glacier r

% Add Regression and Kriged Residuals
sweRK(9).G4 = 9999;
for r = 2:9
    for g = 1:3
        glacier = char(options.glacier(g));
        sweRK(r).(glacier) = sweBMS(r).(glacier) + res_bmaKRIG(r).(glacier).pred;
    end   
end
    
 %Residuals of RK
 sampledRK(9).G4 = 9999;    sampledKRIG(9).G4 = 9999;   sampledBMA(9).G4 = 9999;
 residualsRK(9).G4 = 9999;  residualsKRIG(9).G4 = 9999;
for g = 1:3
    glacier = char(options.glacier(g)); 
T = sub2ind(size(sweRK(2).(glacier)),floor(options.N.(glacier)),floor(options.E.(glacier)));
for r = 2:9
sampledRK(r).(glacier)      = sweRK(r).(glacier)(T); 
sampledKRIG(r).(glacier)    = sweKRIG(r).(glacier).pred(T);
sampledBMA(r).(glacier)     = sweBMS(r).(glacier)(T);

residualsKRIG(r).(glacier)  = SWE(g).swe - sampledKRIG(r).(glacier);
residualsRK(r).(glacier)    = SWE(g).swe- sampledRK(r).(glacier);

sweRK(r).LOO.(glacier)         = res_bmaKRIG(r).LOO.(glacier) + sampledBMA(r).(glacier);
end
end
    clear E g glacier N r T

% clf
% surf(sweRK(2).G13, 'FaceAlpha',0.5,'LineStyle','none'); hold on
% plot(E,N,'.k')

%% Plotting -> residuals
    param = 'pred';
    topoParam.G4  = res_bmaKRIG(4).G4.(param);
    topoParam.G2  = res_bmaKRIG(4).G2.(param);
    topoParam.G13 = res_bmaKRIG(4).G13.(param);
    topoParam.rig = options.rig;

    PlotTopoParameter(topoParam,'symmetric', 'Residual (m w.e.)', SWE, 'symmetric','massB')
        saveFIG('residualsKriged',18,'3G')

%% Plotting -> regression kriging
%     param = 'RK';
%     opt = 8;
%     topoParam.G4  = sweRK(opt+1).G4;
%     topoParam.G2  = sweRK(opt+1).G2;
%     topoParam.G13 = [sweRK(opt+1).G13; nan(10,size(sweRK(opt+1).G13,2))];
%     topoParam.G13 = [sweRK(opt+1).G13, nan(size(sweRK(opt+1).G13,1),10)];
% 
%     PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWE, 'black', 'massB')

    load Full.mat fullRK
    load TopoSWE.mat SWE
PlotTopoParameter(fullRK.S2,'RK', 'SWE (m w.e.)', SWE, 'black', 'massB')   
    saveFIG('RegressionKriging',18,'3G')

    
%% Model param table
nugget = table(zeros(8,1), zeros(8,1),zeros(8,1),...
                'RowNames', options.densityName, 'VariableNames', options.glacier);
maxLL = nugget;  
model = sweKRIG;

for r = 2:9
    for g = 1:3
            glacier = char(options.glacier(g)); 
nugget{r-1,g}   = round(model(r).Model.(glacier).nugget,3);
maxLL{r-1,g}    = round(model(r).Model.(glacier).maxLL,1);   
    end
end
    clear model r g glacier
writetable(nugget, '/Users/Alexandra/Downloads/modelparam.csv','FileType','text')

%% Plot - Actual vs fitted data  
    load TopoSWE.mat topo_sampled 
    load Full.mat fullSWE fullRK options
RKinput = SampledCell(fullRK.S2);
WBinput = ObsInCell(fullSWE.S2.input, topo_sampled);

close all
    dim = [0.16 0.5 0.3 0.3;...
           0.44 0.5 0.3 0.3;...
           0.72 0.5 0.3 0.3];
for i = 1:3;    glacier = options.glacier{i};
    yObserved   = WBinput.(glacier)(:,1);
    yModel      = RKinput.(glacier);
    
    subplot(1,3,i)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
        plot(yObserved, yModel, '.', 'Color', options.RGB(i,:),'MarkerSize',13); hold on

        [f.(glacier), g.(glacier)] = fit(yObserved, yModel,'poly1');
        p = plot(f.(glacier)); hold on
        set(p,'Color',options.RGB(i,:)); set(p, 'LineWidth',1.5);     
        xlabel('Measured WB (m w.e.)'); ylabel('Modelled WB (m w.e.)');
        title(['Glacier ',glacier(2:end)])
                axis square;    box on
        b = gca; legend(b,'off');
        annotation('textbox',dim(i,:),'String', ['R^2=',num2str(round(g.(glacier).rsquare,2), '%.2f')],'FitBoxToText','on')
end

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',12)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 12 4];
 saveFIG('krigRKfit')

    clear b coeffs dim f fig filename g glacier i line option p r X yModel yObserved

%% Plots - Actual vs fitted data for all SWE options
clf
R2value = [];

figure(1)
for i = 1:3

    for r = 2:9
    option = r;
    name	= char(options.glacier(i));
    y       = SWE(i).swe;
    X       = sweRK(r).LOO.(name);
   
    subplot(1,3,i)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
        [f.(name), g.(name)] = fit(y, X,'poly1');
        p = plot(f.(name)); hold on
        set(p,'Color',options.RGB(i,:)); set(p, 'LineWidth',1.5);     
        xlabel('Measured Winter Balance (m w.e.)'); ylabel('Modelled Winter Balance (m w.e.)');
        title(name)
                axis square;    box on        
        R2value = mean([R2value g.(name).rsquare]);
        if r == 8; p.Color = [0 0 0]; end
    end
        b = gca; legend(b,'off');
        dim = [b.Position(1)+0.01 b.Position(2)+.37 .3 .3];
        annotation('textbox',dim,'String', ['R^2=',num2str(round(R2value,2), '%.2f')],'FitBoxToText','on')
    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 12 4];

end
saveFIG('sweRKfit_allLines')
        clear b p R2value dim f g j r fig filename h i params X y

%% Plotting -> residual distribution (kriged) as % of regressed SWE    

opt = 8;
for g = 1:3
    glacier = char(options.glacier(g));
residualPer.(glacier) = res_bma_KRIG(opt).(glacier).pred./sweBMS(opt).(glacier)*100;

residualPer.(glacier)(isinf(residualPer.(glacier))) = NaN;
residualPer.(glacier)(residualPer.(glacier) < -100) = -100;
residualPer.(glacier)(residualPer.(glacier) > 100)  = 100;

residualPer.(glacier) = abs(residualPer.(glacier));
end
        
    PlotTopoParameter(residualPer, 'residual%', {'Kriged residuals as percent', 'of BMA estimate (%)'},...
                        SWE, 'black', 'NOmassB')
        saveFIG('Residuals_Percent_BMA','3G')

%% Plotting -> all opts kriged SWE map

glacier = 'G2';
for i = 1:length(options.densityName)
subplot(4,2,i)
imagesc(sweKRIG(i+1).(glacier).pred)
    caxis([0 1.5])
    axis equal
    title([options.densityName(i), num2str(nanmean(sweKRIG(i+1).(glacier).pred(:)))])
end


%% Plot -> CI for kriging
load Full.mat fullSK options
load TopoSWE.mat sweOPT

opt = 'S2';
for g = 1:3;   glacier = options.glacier{g};
   CI.(glacier) = (fullSK.(opt).(glacier).upper95 - fullSK.(opt).(glacier).lower95)...
       ./fullSK.(opt).(glacier).pred*100;
   
   CI.(glacier)(isinf(CI.(glacier))) = NaN;
   %CI.(glacier)(CI.(glacier)>400) = 400;
end

PlotTopoParameter(CI,'uncertainity', 'Confidence interval (%)', ...
                    sweOPT(2), 'black', 'NmassB')
    saveFIG('KrigingCI_percent',18,'3G')
    
    
    
data = abs(SWE7(2).swe - SWE9(2).swe);
TT = data>0.01;
scatter(SWE(2).utm(TT,1), SWE(2).utm(TT,2),20, data(TT),'filled'); colorbar
TT = data>0.02;
scatter(SWE(2).utm(TT,1), SWE(2).utm(TT,2),20, data(TT),'filled'); colorbar
figure
plot(SWE7(2).swe(TT),'.')
hold on
plot(SWE9(2).swe(TT),'.')
SWE(2).label(TT)


%% Cross correlation with elevation

[acor,lag] = xcorr(SWE(3).swe,topo_sampled_ns.G13.elevation);

[~,I] = max(abs(acor));
lagDiff = lag(I)
    
%% Variograms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear
close all
load TopoBMS_MLR.mat

for g = 1:3
    glacier = char(options.glacier(g));

x = SWE(g).utm(:,1);%-min(SWE(g).utm(:,1));
y = SWE(g).utm(:,2);%-min(SWE(g).utm(:,2));
z = residualsBMS(8).(glacier);
%z = SWE(g).swe;

% Create data pairs
% Construct arrays that include all combinations of data pairs
  [X1,X2] = meshgrid(x);
  [Y1,Y2] = meshgrid(y);
  [Z1,Z2] = meshgrid(z);
  
% Compute separation distances between all pairs:
  D = sqrt( (X1-X2).^2 + (Y1-Y2).^2 );
  
% Compute semivariance for each data pair (this counts each pair twice):
  G = (1/2)*(Z1 - Z2).^2;
  
% Define minimum and maximum lags for irregularly spaced data

% Min lag: mean min distance of points
% Assign diagonal of D to value of NaN (to remove zeros)    
  D2 = D + D.*diag(x*NaN);  
  
  lag = mean(min(D2));
  
% Max lag: half max distance between points
  hmd = max(max(D))/2;           % Find max distance between points/2
  max_lags = floor(hmd/lag);     % Calculate a number of lags
  
% Bin calculated distances as integers according to lag  
  all_lags = ceil(D/lag);
  
  for i=1:max_lags
     selection = (all_lags == i);               % select points in bin (1 = true)
     mean_lag(i) = mean(mean(D(selection)));    % compute mean lag
     num_pairs(i) = sum(sum(selection == 1))/2; % count pairs
     vario(i) = mean(mean(G(selection)));       % populate variogram
  end
  
% Compute variance of original data
  varz = var(z);

% Plot semivariogram and variance of data 
% figure(1)
% subplot(1,3,g)
% plot(mean_lag,vario,'o','Markerfacecolor','b')
%   hold on
%   plot([0 max(mean_lag)],[varz varz],'--k')
%   grid on
%   ylim([0 1.05*max(vario)]);
%   xlabel('Lag')
%   ylabel('Semivariance')
%   title(glacier)
%       saveFIG('Variogram_BMSresiduals');



% Kriging**
    glacier = char(options.glacier(g));  
% Create lags vector for model that extends to max expt lag
  lags = [0:1:ceil(max(mean_lag))];
  
% Spherical model with optional nugget:
  nugget = min(vario);
  sill = varz;
  range = 500;
  
%   mod_sphere = nugget + ...
%                ( sill*( 1.5*(lags/range) - 0.5*(lags/range).^3).*...
%                (lags <= range) + ...
%                sill*(lags>range));
%   hold on
%   plot(lags,mod_sphere,'k','Linewidth',1.25)
% 
% %% Linear model with nugget:
%  
%   slope = (0.007-0.0025)/500;
%   
%   mod_lin = (nugget + slope*lags).*(nugget + slope*lags < varz) + ...
%                      varz*(nugget + slope*lags >= varz);
% 
%   hold on
%   plot(lags,mod_lin,'g','Linewidth',1.25)
% 
% %% Exponential model with nugget:
%  
%   mod_exp = nugget + sill.*(1 - exp(-3*lags/range));
% 
%   hold on
%   plot(lags,mod_exp,'c','Linewidth',1.25)
%   legend('data','variance','spherical model','linear model',...
%          'exponential model','Location','SouthEast')

% Populate semivariance matrix W based on exponential model
  W = nugget + sill.*(1 - exp(-3*D/range));

% Add extra row and column to W for Lagrange multiplier
  N = length(x);        
  W(:,end+1) = 1;     % Add final column of ones
  W(end+1,:) = 1;     % Add final row of ones
  W(end,end) = 0; % Add corner value of zero
  
% Invert W matrix for later
  Winv = inv(W);
  
% Construct regular grid over which interpolated values are required
  H = size(topo_full.(glacier).elevation);
  gridS = 40;
  xgrid = [min(rig.(glacier)(:,1)):gridS:max(rig.(glacier)(:,1))];
  ygrid = [min(rig.(glacier)(:,2)):gridS:max(rig.(glacier)(:,2))];
  
  [Xgrid,Ygrid] = meshgrid(xgrid,ygrid);
  
% Convert arrays to single vectors:
  Xvec = reshape(Xgrid,[],1);
  Yvec = reshape(Ygrid,[],1);
  
% Initialize variables
  Zvec = 99999*ones(1,length(Xvec));
  S2vec = 99999*ones(1,length(Xvec));
  
% Computed kriged estimate at each point i
  for i=1:length(Xvec)
    % calculate distances between point of interest and all observations:
      D0 = sqrt( (x - Xvec(i)).^2 + (y - Yvec(i)).^2 ); 
    % calculate B (RHS) based on distances above and variogram model
      B = nugget + sill.*(1 - exp(-3*D0/range));
    % Add final element for Lagrange multiplier
      B(N+1) = 1;
    % Compute kriging weights and lagrange multiplier with matrix multiplication
      Lambda = Winv*B;
    % Compute kriging estimate at point of interest as weighted sum of obs    
      Zvec(i) = z'*Lambda(1:N);    
            % Zvec(i) = sum(Lambda(1:N).*z);
    % Compute kriging variance at point of interest
      S2vec(i) = B'*Lambda; 
            % S2vec(i) = sum(Lambda(1:N).*B(1:N)) + Lambda(N+1);
  end
  
% Reshape results and plot  
  Mx  = length(xgrid); My = length(ygrid);
  Z  = reshape(Zvec,My,Mx); 
        dataN = flipud(topo_full.(glacier).elevation);
        sizeZ = size(Z);    sizedataN = size(dataN);
        if sizeZ(1,1) ~= sizedataN(1,1)
           Z = Z(2:end,:); 
        elseif sizeZ(1,2) ~= sizedataN(1,2)
           Z = Z(:,2:end); 
        end
        Z(isnan(dataN)) = NaN;
  S2 = reshape(S2vec,My,Mx);

figure(2)
subplot(1,3,g)
  imagesc(xgrid,ygrid,Z);  hold on
  plot(x,y,'w.'); hold on
  plot(rig.(glacier)(:,1),rig.(glacier)(:,2),'k')
  axis('xy')
  axis image
  xlabel('x coordinate')
  ylabel('y coordinate') 
  colorbar('vert')
  title('Kriged estimate z_{hat}')
  end
%   figure
%   imagesc(xgrid,ygrid,S2)
%   hold on
%   plot(x,y,'w.')
%   axis('xy')
%   axis image
%   xlabel('x coordinate')
%   ylabel('y coordinate') 
%   colorbar('vert')
%   title('Kriging variance with observation points')
% **********************************************************
  

% figure(1)
%     saveFIG('Variogram_BMSresiduals','3G');
% figure(2)
%     saveFIG('Kriging_BMSresiduals');






%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Kriging 1000 reps with 2/3 data

clear; close all
load TopoBMS_MLR.mat sweOPT options

nRuns = 100;
    KRIGsubset(nRuns+1).G4 = 9999; 
    KRIGmodel(nRuns+1).G4 = 9999; 
    outD(nRuns+1).G4 = 9999; 

for d = 1:8; den = options.DenOpt{d};    
for i = 1:nRuns
 for g = 1:3;     glacier = options.glacier{g};

    fullD = sweOPT(d+1).(glacier);
        
    %Random subsets (2/3)
    I.(glacier) = randi(length(fullD),floor(length(fullD)*2/3),1);
    Inot.(glacier) = ~ismember(1:length(fullD),I.(glacier)); %%%

    inputD = fullD(I.(glacier),:);

    %Kriging
    outD(i).(den).(glacier) = KrigingR(inputD(:,1), inputD(:,2:3), glacier);
    
    KRIGsubset(i).(den).(glacier) = outD(i).(den).(glacier).pred;
    KRIGmodel(i).(den).(glacier) = outD(i).(den).(glacier).Model.theta;  
 end    
     estD = SampledCell(KRIGsubset(i).(den));
   
    %RMSE  
   for g = 1:3;     glacier = options.glacier{g};
       obsD = sweOPT(d+1).(glacier)(Inot.(glacier),1);
       predD = estD.(glacier)(Inot.(glacier));
     RMSE(i).(den).(glacier) = sqrt(mean((predD-obsD).^2));
   end
end

%Get min RMSE option for each glacier
    for i = 1:length(RMSE)
    temp = struct2cell(RMSE(i).(den)(:))';
    rmseD(i,1:3) = temp{:,:};
    end

   for g = 1:3;     glacier = options.glacier{g};
    rmseI = rmseD{:,g}== min(rmseD{:,g});
    WBSK.(den).(glacier) = outD(rmseI).(den).(glacier);
   end


%STD from upper and lower 95%
    for g = 1:3;     glacier = options.glacier{g};
    WBmin = nanmean(WBSK.(den).(glacier).lower95(:));
    WBmax = nanmean(WBSK.(den).(glacier).upper95(:));
    Qstd.(den).(glacier) = (WBmax-WBmin)/2/1.96;
    end

end
  
  
  
%% Universal Kriging

%clear; load LR_SK_RK.mat sweSK sweUK 
OPTIONS
%SK = sweSK(2:9);  UK = sweUK(2:9);    clear swe*

%Plot diff between UK and SK (UK minus SK)
sd = [1 5 2 6 3 7 4 8];
Cmap = flipud(cbrewer('div','RdBu',20,'PCHIP'));
for g = 1:3; glacier = options.glacier{g};
   figure(3);
    for d = 1:8; D = options.DenOpt{d};
    subplot(2,4,sd(d))
        Pdata = UK(d).(glacier).pred - SK(d).(glacier).pred;
        h = imagesc(Pdata); colorbar
            set(h,'alphadata',~isnan(Pdata));    axis off
            title(D)
            colormap(Cmap)
    end
    saveFIG(['UKminusSK_',glacier])
end
  
%Plot UK WB estimate
sd = [1 5 2 6 3 7 4 8];
for g = 1:3; glacier = options.glacier{g};
   figure(3);
    for d = 1:8; D = options.DenOpt{d};
    subplot(2,4,sd(d))
        Pdata   = UK(d).(glacier).pred;
        h = imagesc(Pdata); colorbar
            set(h,'alphadata',~isnan(Pdata));    axis off
            title([D,' (',num2str(round(nanmean(nanmean(Pdata)),2)),' m w.e.)'])
    end
    saveFIG(['UKestimatedWB_',glacier])
end
  
 %% 
%Derek Plots
GDen = [6,5];% corresponds to ['F3';'S3'] on glaciers 2 and 13;
Cmap_diff   = flipud(cbrewer('div','RdBu',100,'PCHIP'));
Cmap_wb     = cbrewer('seq','Greys',100,'PCHIP');
   figure(3); clf; n=1;
for g = 2:3; glacier = options.glacier{g};
    d = GDen(g-1);
    s1 = subplot(2,3,n);
        Pdata   = SK(d).(glacier).pred;
        h = imagesc(Pdata); Cbar = colorbar;
            set(h,'alphadata',~isnan(Pdata));    axis off
            title(['SK ',glacier,' (',num2str(round(nanmean(nanmean(Pdata)),2)),' m w.e.)'])
            colormap(s1, Cmap_wb);     caxis([0 0.55]); 
            ylabel(Cbar,'Winter balance (m w.e.)')        

    s2 = subplot(2,3,n+1);
        Pdata   = UK(d).(glacier).pred;
        h = imagesc(Pdata); Cbar = colorbar;
            set(h,'alphadata',~isnan(Pdata));    axis off
            title(['OK ',glacier,' (',num2str(round(nanmean(nanmean(Pdata)),2)),' m w.e.)'])
            colormap(s2, Cmap_wb);    caxis([0 0.55])        
            ylabel(Cbar,'Winter balance (m w.e.)')        

    s3 = subplot(2,3,n+2);
        Pdata   = UK(d).(glacier).pred-SK(d).(glacier).pred;
        h = imagesc(Pdata); Cbar = colorbar;
            set(h,'alphadata',~isnan(Pdata));    axis off
            title(['OK minus SK ',glacier])
            colormap(s3, Cmap_diff)
            ylabel(Cbar,'Winter balance (m w.e.)')        
    n=4;        
end
    saveFIG('SKandOK_G2_G13')
