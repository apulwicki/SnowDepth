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
        text(SP(g,1)-1, FS(g,1)-4, strG,...
            'HorizontalAlignment','right','FontSize',9,'FontName','Arial');
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

%clear; load WSMBDistribution.mat
figure(9); clf;
x = 0.25:0.001:0.75;
for o = 1:3
    if      o == 1; data = varBzz;     Pmax = 53;   t = '\sigma_{ZZ} Variability';       f = 'zz';
    elseif  o == 2; data = varBbeta;   Pmax = 20.3;   t = '\sigma_{\beta} Variability';             f = 'beta';
    elseif  o == 3; data = varBbetazz; Pmax = 18.8;   t = '\sigma_{\beta} and \sigma_{ZZ} Variability'; f = 'betaNzz';
    end
for d = 1:8;    den = options.DenOpt{d};
for g = 1:3;    glacier = options.glacier{g};
    
ProbDen.(den).(glacier) = fitdist(data.(den).(glacier)(:),'Normal');
    y = pdf(ProbDen.(den).(glacier),x);   y = y/Pmax;  

subplot(1,3,o) 
fill(x,y,options.RGB(g,:),'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)');  
    
end
end
    if o ==1; legend(options.glacier,'location','northwest')
    end
    xlim([min(x) max(x)])
    title(t);
end

saveFIG_IGS('WSMBDist_LR',2,8.6)

%% WSMB Distribution - LR vs SK

%clear; load WSMBDistribution.mat
figure(10); clf;
x = 0.2:0.001:0.75;

C = [115 191 184; 255 144 0]/255;
G = {'Glacier 4','Glacier 2','Glacier 13'};

Pmax = [76, 18, 37];

for g = 1:3;    glacier = options.glacier{g};  
for d = 1:8;    den = options.DenOpt{d};
    
ProbDenLR.(den).(glacier) = fitdist(varBbetazz.(den).(glacier)(:),'Normal');
    yLR = pdf(ProbDenLR.(den).(glacier),x);   yLR = yLR/Pmax(g);  
ProbDenSK.(den).(glacier) = fitdist(varBkriging.(den).(glacier)(:),'Normal');
    ySK = pdf(ProbDenSK.(den).(glacier),x);   ySK = ySK/Pmax(g);  
    
subplot(1,3,g) 
fill(x,yLR,C(1,:),'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
fill(x,ySK,C(2,:),'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
    
end

    if g ==1; legend([{'LR'},{'SK'}],'location','northwest')
    end
    title(G{g});
end

saveFIG_IGS('WSMBDist_LRvsSK',2,8.6)

%% WSMB Distribution -> full PDF LR and SK

%clear; load WSMBDistribution.mat varBbetazz varBkriging

for g = 1:3
glacier = options.glacier{g};
for d = 1:8
den = options.DenOpt{d};
Tdbetazz.(glacier)(:,:,d)       = varBbetazz.(den).(glacier);
Tdzzkriging.(glacier)(:,:,d)    = varBkriging.(den).(glacier);
end
end

figure(12); clf
x = 0.2:0.001:0.75;
Pmax = [13.8, 24];
for g = 1:3;     glacier = options.glacier{g};

ProbDenLR.(glacier) = fitdist(Tdbetazz.(glacier)(:),'Normal');
    yLR = pdf(ProbDenLR.(glacier),x);   yLR = yLR/Pmax(1);  
ProbDenSK.(glacier) = fitdist(Tdzzkriging.(glacier)(:),'Normal');
    ySK = pdf(ProbDenSK.(glacier),x);   ySK = ySK/Pmax(2);  
    
subplot(2,1,1) 
fill(x,yLR,options.RGB(g,:),'FaceAlpha',0.7, 'EdgeColor', 'none'); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
subplot(2,1,2)     
fill(x,ySK,options.RGB(g,:),'FaceAlpha',0.7, 'EdgeColor', 'none'); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
   
end
ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
legend(options.glacier,'Location','best')

saveFIG_IGS('WSMBDist_full',1,13)

%% WSMB Distribution - spatial variability (one density)

%clear; load WSMBDistribution.mat D
load TopoSWE.mat SWE

% Linear regression
den = 'S2';
for g = 1:3; glacier = options.glacier{g};
    DD.(glacier) = D_LR.(den).(glacier)/(max(D_LR.(den).(glacier)(:))*0.55);
    DD.(glacier)(options.mapNaN.(glacier)) = NaN;
end

figure(11); 
PlotTopoParameter_IGS(DD,'hot','Variability',SWE,'none','nomassB')
    saveFIG_IGS('SpatialVar_LR',2,8.6)

% Simple krigings
den = 'S2';
for g = 1:3; glacier = options.glacier{g};
    DD.(glacier) = D_SK.(den).(glacier)/(max(D_SK.(den).(glacier)(:))*0.9);
    DD.(glacier)(options.mapNaN.(glacier)) = NaN;
end

figure(12); 
PlotTopoParameter_IGS(DD,'hot','Variability',SWE,'none','nomassB')
    saveFIG_IGS('SpatialVar_SK',2,8.6)

