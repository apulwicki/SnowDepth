%% Measurement - Depth Boxplot
    clear
load TopoSWE.mat allDepth options

    lG4     = length(allDepth.G4(:)); 
    lG2     = length(allDepth.G2(:)); 
    lG13    = length(allDepth.G13(:));
    toBox = [[allDepth.G4(:); nan(lG13-lG4,1)],...
             [allDepth.G2(:); nan(lG13-lG2,1)],...
              allDepth.G13(:)];
figure(3); clf
boxplot(toBox,'labels',{'Glacier 4','Glacier 2','Glacier 13'},...
               'BoxStyle','outline','Colors',options.RGB,'Width',0.35,...
               'OutlierSize',4,'Symbol','o')
    ylabel('Snow depth (cm)');
saveFIG_IGS('DepthBoxplot',1,8.6)

%% Density Interp - SP vs FS
    clear
load TopoSWE.mat Density

    den     = [nan(6,1), cell2mat(Density.pitANDtube(:,2:10))];
    SP      = den(:,7);
    FS      = den(:,2);
    errorSP = [den(:,7)-den(:,9),den(:,10)-den(:,7)];   %min and max SP
    errorFS = [den(:,2)-den(:,4),den(:,5)-den(:,2)];    %min and max FS

figure(4); clf
errorbarxy(SP,FS,errorSP(:,2),errorFS(:,2),errorSP(:,1),errorFS(:,1),...
                'Color','k','LineStyle','none','Marker','s',...
                'MarkerFaceColor','k','LineWidth',1,'MarkerSize',7); hold on
    axis([220 400 220 400])
    grid on
    line = refline(1,0);
        line.Color = 'k'; line.LineStyle = '--'; hold on
    xlabel('Snow pit-derived integrated density (kg m^{-3})')
    ylabel('Federal Sampler-derived density (kg m^{-3})')
     
    %Label points
    labels = {'G4\_USP';'G4\_LSP';'G2\_USP';'G2\_LSP';'G13\_ASP';'G13\_USP'};
    for g = 1:length(SP)
        strG = labels{g,1};
        if g ==1
        text(SP(g,1)+45, FS(g,1)-6, strG,...
            'HorizontalAlignment','right','FontSize',9,'FontName','Arial');            
        else
        text(SP(g,1)-3, FS(g,1)+7, strG,...
            'HorizontalAlignment','right','FontSize',9,'FontName','Arial');
        end
    end
    ax = gca; ax.XTick = 220:40:400; ax.YTick = 220:40:400;

saveFIG_IGS('SPvsFS',1,8.6);

%% Grid Cell - Zigzag histrogram
    clear
load TopoSWE.mat SWEzz
    
%     labels(:,1) = {'G4 LZ';'G4 MZ';'G4 UZ';'none'}; 
%     labels(:,2) = {'G2 LZ';'G2 MZ';'G2 UZ';'none'};
%     labels(:,3) = {'G13 LZ';'G13 M_1Z';'G13 M_2Z';'G13 UZ'};
    labels(:,1) = {'L';'M';'U';'none'}; 
    labels(:,2) = {'L';'M';'U';'none'};
    labels(:,3) = {'L';'M_1';'M_2';'U'};

figure(5); clf
    for g = 1:3
        zz = categories(SWEzz(g).ZZ);

    subplot(3,1,g)
        for j = 1:length(zz)
            ZZdata = SWEzz(g).swe(SWEzz(g).ZZ==char(zz(j)));
            ZZdata = (ZZdata-mean(ZZdata));%/mean(ZZdata)*100;%/std(ZZdata);

                if j ==1;   bins    = 20;   %round(sqrt(length(SWEzz(g).swe)));
                            edges   = linspace(-0.15,0.15,bins); end
                            N       = histcounts(ZZdata,edges);   
            plot((edges(:,1:end-1)+edges(:,2:end))/2,N/sum(N),'LineWidth',2); hold on 
            xlabel('Zigzag SWE distribution (m w.e.)');     ylabel('Probability')
            grid on
            xlim([-0.15 0.15])
            %title(options.glacier{g})
            %ax = gca; ax.XTick = [-40:20:40];
        end
            legend(labels{1:length(zz),g},'Location','northeast')
    end
saveFIG_IGS('ZigzagHistogram',1,17.8);

%% Interp Method - LR & SK map
    clear
load Full.mat fullLR fullSK
load TopoSWE.mat SWE

figure(6); clf
PlotTopoParameter_IGS(fullLR.S2, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black', 'massB')
	saveFIG_IGS('LR_map',2,8.6)
figure(6); clf
PlotTopoParameter_IGS(fullSK.S2, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black', 'massB')
	saveFIG_IGS('SK_map',2,8.6)

%% Interp Method - Observed vs Estimated SWE
    clear
load Full.mat fullLR fullSK
load TopoSWE.mat SWE topo_sampled options

den = 'S2';
    yObserved   = ObsInCell(SWE, topo_sampled);

figure(7); clf
for y = 1:2
    if      y ==1;  yEstimated = SampledCell(fullLR.(den));   k = 0;
    elseif  y ==2;  yEstimated = SampledCell(fullSK.(den));   k = 3;
    end

for g = 1:3;    glacier = options.glacier{g};
    
subplot(2,3,g+k)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
    plot(yObserved(g).swe, yEstimated.(glacier), 'o', 'Color', options.RGB(g,:),'MarkerSize',2); hold on

        [F.(glacier), G.(glacier)] = fit(yObserved(g).swe, yEstimated.(glacier),'poly1');
        p = plot(F.(glacier)); hold on
            set(p,'Color',options.RGB(g,:)); set(p, 'LineWidth',1.5);     
        xlabel('Observed (m w.e.)'); 
        if      y ==1;  ylabel('LR Estimate (m w.e.)');
        elseif  y ==2;  ylabel('SK Estimate (m w.e.)');
        end
        %title(['Glacier ',glacier(2:end)])
                axis square;    box on;     grid on
        b = gca; legend(b,'off');
        %dim = [b.Position(1)+0.01 b.Position(2)+.37 .3 .3];
        %annotation('textbox',dim,'String', ['R^2=',num2str(round(G.(glacier).rsquare,2))],'FitBoxToText','on')
end
end

saveFIG_IGS(['observedVSestimated_',den],2,10)

%% Interp Method - Beta coeffs boxplot
    
    clear
load Full.mat fullLR options
load TestInterp.mat SigBetaDen

%Get betas for all density options
betas = zeros(8,7,3);
for d = 1:8; den = options.DenOpt{d};
    betas(d,:,:) = fullLR.(den).coeff{1:7,:}; 
end

for g = 1:3;    glacier = options.glacier{g};
    betas(:,:,g) = betas(:,:,g).*SigBetaDen.(glacier)(2:end,:)';
end

BETAS = [reshape(betas(:,:,1), [1 size(betas(:,:,1))]);...
         reshape(betas(:,:,2), [1 size(betas(:,:,2))]);...
         reshape(betas(:,:,3), [1 size(betas(:,:,3))])];

figure(8); clf        
aboxplot(BETAS,'labels',options.topoVars, ...
         'Colormap', options.RGB, 'OutlierMarkerSize',6); % Advanced box plot
        legend('Glacier 4','Glacier 2','Glacier 13'); % Add a legend
        ylabel('Regression coefficient'); hold on 
        
        ylim([-0.06 0.12])
        line(xlim,[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        for i = 1:6
        line([i+0.5 i+0.5],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end
        
saveFIG_IGS('BetaCoeffs',1,8.6)

%% WSMB Distribution - LR sources of var

%clear; load WSMBDistribution.mat var* options

%get std of kriging surfaces
for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};
    
    for i = 1:length(tempSKzz.(den))
Qmin(i) = nanmean(tempSKzz.(den)(i).(glacier).lower95(:));
Qmax(i) = nanmean(tempSKzz.(den)(i).(glacier).upper95(:));
Qpre(i) = nanmean(tempSKzz.(den)(i).(glacier).pred(:));
    end
Qstd.(den).(glacier) = (Qmax-Qmin)/2/1.96;

Qstd.(den).(glacier) = mean(Qstd.(den).(glacier));
Qmean.(den).(glacier) = mean(Qpre);

end
end
    

figure(9); clf;
x = 0:0.001:1.1;
for o = 1:4
    if      o == 1; data = varBzz;     Pmax = 53;   t = '\sigma_{SWE} Variability';     f = 'zzLR';
    elseif  o == 2; data = varBbeta;   Pmax = 20.3;   t = '\sigma_{\beta} Variability'; f = 'beta';
    elseif  o == 3; data = varBkriging; Pmax = 76;   t = '\sigma_{SWE} Variability';  f = 'zzSK';
    elseif  o == 4; data = varBkriging; Pmax = 4.1;   t = '\sigma_{KRIG} Variability'; f = 'KRIG';
    end
for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};
    
if      o <= 3    
ProbDen.(den).(glacier) = fitdist(data.(den).(glacier)(:),'Normal');
    y = pdf(ProbDen.(den).(glacier),x);   y = y/Pmax;
elseif  o == 4
ProbDen.(den).(glacier) = makedist('Normal',Qmean.(den).(glacier),Qstd.(den).(glacier));
    y = pdf(ProbDen.(den).(glacier),x);   y = y/Pmax;
    y = [0 y]; x = [x 0];
end


subplot(2,2,o) 
p(g) = fill(x,y,options.RGB(g,:),'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)');  
    
end
end
    if o ==2; legend(p,options.glacier,'location','northeast')
    end
    xlim([min(x) max(x)])
    title(t);
end

saveFIG_IGS('WSMBDist_LR',2,12)


%% WSMB Distribution -> full PDF LR and SK

%clear; load WSMBDistribution.mat varBbetazz varBkriging Q* options
for g = 1:3
glacier = options.glacier{g};
for d = 1:8
den = options.DenOpt{d};
Tdbetazz.(glacier)(:,:,d)       = varBbetazz.(den).(glacier);
Tdzzkriging.(glacier)(:,:,d)    = varBkriging.(den).(glacier);
end
end

%get std of kriging surfaces
for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};
    
    for i = 1:length(tempSKzz.(den))
Qmin(i) = nanmean(tempSKzz.(den)(i).(glacier).lower95(:));
Qmax(i) = nanmean(tempSKzz.(den)(i).(glacier).upper95(:));
Qpre(i) = nanmean(tempSKzz.(den)(i).(glacier).pred(:));
    end
Qstd.(den).(glacier) = (Qmax-Qmin)/2/1.96;
Qstd.(den).(glacier) = mean(Qstd.(den).(glacier));
Qmean.(den).(glacier) = mean(Qpre);

Qstd.(den).(glacier) = mean(Qstd.(den).(glacier));
Qmean.(den).(glacier) = mean(Qpre);

Tsk.(glacier)(:,:,d) = normrnd(Qmean.(den).(glacier),Qstd.(den).(glacier),100);
end
end

figure(12); clf
x = 0:0.001:1.1;
Pmax = [13.8, 24.5];
for g = 1:3;     glacier = options.glacier{g};

ProbDenLR.(glacier) = fitdist(Tdbetazz.(glacier)(:),'Normal');
    yLR     = pdf(ProbDenLR.(glacier),x);   %yLR = yLR/Pmax(1);  
    YlrN    = yLR./trapz(x,yLR);
ProbDenSK.(glacier) = makedist('Normal',T(g,3),T(g,4));
    ySK     = pdf(ProbDenSK.(glacier),x);   %ySK = ySK/Pmax(2);
    YskN    = ySK./trapz(x,ySK);
    
subplot(2,1,1) 
fill(x,yLR,options.RGB(g,:),'FaceAlpha',0.85, 'EdgeColor', 'none'); hold on
    ylabel('LR Probability'); xlabel('WSMB (m w.e.)'); 
    legend(options.glacier,'Location','best')
subplot(2,1,2)     
fill([0 x],[0 ySK],options.RGB(g,:),'FaceAlpha',0.85, 'EdgeColor', 'none'); hold on
    ylabel('SK Probability'); xlabel('WSMB (m w.e.)'); 

    xlim([min(x),max(x)])
end

saveFIG_IGS('WSMBDist_full',1,10)

%% WSMB Distribution - total spatial variability 

%clear; load WSMBDistribution.mat D*
load TopoSWE.mat SWE

for o = 2%1:2
    if o == 1
    D = D_LR;   s = 'LR';   m = 0.6;
    elseif o == 2
    D = D_SK;   s = 'SK';   m = 0.95;
    end
        

% Linear regression
for g = 1:3; glacier = options.glacier{g};
    DPlot.(glacier) = zeros(options.mapsize(g,:));
    for d = 1:8; den = options.DenOpt{d};
        DD.(glacier) = D.(den).(glacier)/(max(D.(den).(glacier)(:)));
        DD.(glacier)(options.mapNaN.(glacier)) = NaN;
    DPlot.(glacier) = DPlot.(glacier) + DD.(glacier);
    end
    DPlot.(glacier) = DPlot.(glacier)/nanmax(DPlot.(glacier)(:)*m);
end

figure(o);
PlotTopoParameter_IGS(DPlot,'hot','Variability',SWE,'none','nomassB')
    saveFIG_IGS(['SpatialVar_',s],2,8.6)
end

%% Accumulation gradient

taylor(:,1) = [571731.48;577258.98;580978.1;587346.4;591126.5;597353.2;601796.1;608101];
taylor(:,2) = [6737517.35;6733918.68;6730286.9;6730436.4;6724959.2;6730694.1;6734532.2;6736574.4];
taylor(:,3) = [2620;2640;2380;2225;2070;1915;1765;1615];
taylor(:,4) = [3.2;3.7;2.6;2.7;2.4;1.78;1.65;0.91];
taylor(:,5) = [0.407;0.409;0.394;0.364;0.388;0.385;0.390;0.342];
taylor(:,6) = [1.302;1.513;1.024;0.983;0.932;0.685;0.643;0.311];

alex(:,1) = [566453.4;570077.4;595349.9;601424.8;605031];
alex(:,2) = [6727621.2;6732429.7;6741065.2;6753607.6;6762717.2];
alex(:,3) = [2610;2730;2321;2472;2434];
alex(:,6) = [1.30;1.59;0.5844;0.5785;0.3834];
    alexerr = [0.029,0.049,0.032]*1.96;
    
    Dt = sqrt((taylor(:,1)-alex(1,1)).^2+(taylor(:,2)-alex(1,2)).^2)/1000;
    Da = sqrt((alex(:,1)-alex(1,1)).^2+(alex(:,2)-alex(1,2)).^2)/1000;
    [Dat, Iat] = sort([Da; Dt]); AT = [alex(:,6); taylor(:,6)]; AT = AT(Iat);    
[Fat, Gat] = fit(Dat,AT,'poly1');
    CIat = predint(Fat,Dat);    
    
figure(2); clf
    pt = plot(Fat);    set(pt,'Color','k'); set(pt, 'LineWidth',1.5); hold on
%     plot(Dat, CIat(:,1),'--', 'Color','k', 'LineWidth',.2)
%     plot(Dat, CIat(:,2),'--','Color','k', 'LineWidth',.2)
L(2) = plot(Da(1:2),alex(1:2,6),'.', 'MarkerSize',13, 'Color',[234, 140, 46]/255);
    errorbar(Da(3:5),alex(3:5,6),alexerr,'.', 'MarkerSize',13, 'Color',[234, 140, 46]/255)
L(1) = plot(Dt,taylor(:,6),'.', 'MarkerSize',13, 'Color',[68, 181, 226]/255); hold on
    xlim([min(Dat), max(Dat)+1]); 
    ylim([0 2])
    grid on
        xlabel('Distance from mountain divide (km)'); ylabel('SWE (m w.e.)')
        legend(L,{'Taylor-Barge (1969)','Pulwicki et al. (2017)'},'Box','off', 'Location','northoutside')
saveFIG_IGS('AccumGrad',1,8.6)
