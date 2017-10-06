%% Measurement - Depth Boxplot and SP vs FS
    clear; close all
load TopoSWE.mat allDepth Density
run OPTIONS 

 %SP vs FS
    den     = [nan(6,1), cell2mat(Density.pitANDtube(:,2:10))];
    SP      = den(:,7);
    FS      = den(:,2);
    errorSP = [den(:,7)-den(:,9),den(:,10)-den(:,7)];   %min and max SP
    errorFS = [den(:,2)-den(:,4),den(:,5)-den(:,2)];    %min and max FS
    markerC = [options.RGB(1,:);options.RGB(1,:);...
               options.RGB(2,:);options.RGB(2,:);...
               options.RGB(3,:);options.RGB(3,:)]; 
    markerS = {'s','o','s','o','^','s'};

for i = 1:length(SP)
errorbarxy(SP(i),FS(i),errorSP(i,2),errorFS(i,2),errorSP(i,1),errorFS(i,1),...
                'Color','k','LineStyle','none','Marker',markerS{i},...
                'MarkerFaceColor',markerC(i,:),'LineWidth',1,'MarkerSize',8,...
                'MarkerEdgeColor','none'); hold on
end
    axis([220 400 220 400])
    grid on
    line = refline(1,0);
        line.Color = 'k'; line.LineStyle = '--'; hold on
    xlabel('Snow-pit-derived density (kg m^{-3})')
    ylabel('Federal-Sampler-derived density (kg m^{-3})')
     
    %Label points
    labels = {'G4\_Mid';'G4\_Low';'G2\_Mid';'G2\_Low';'G13\_High';'G13\_Mid'};
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
        
saveFIG_IGS('SPvsFS',1,7.5)


%% Interp Method - LR map
    clear
load Full.mat fullLR 
load TopoSWE.mat SWE

figure(6); clf
PlotTopoParameter_IGS(fullLR.S2, 'modelledSWE', 'WB (m w.e.)', SWE, 'black', 'massB')
	saveFIG_IGS('LR_map',2,8.6)


%% Interp Method - LR Observed vs Estimated SWE
    clear
load Full.mat fullLR fullSK
load TopoSWE.mat SWE topo_sampled
run OPTIONS.m

den = 'S2';
    yObserved   = ObsInCell(SWE, topo_sampled);
    for g = 1:3;    glacier = options.glacier{g};
    SKinput.(glacier) = fullSK.(den).(glacier).pred;
    end
    
          locX = [.15 .43 .71]; locX = [locX locX];
          locY = [.83; 0.35];     locY = repmat(locY,1,3); locY = [locY(1,:) locY(2,:)];

figure(7); clf
yEstimated = SampledCell(fullLR.(den)); 
   
for g = 1:3;    glacier = options.glacier{g};
    
subplot(1,3,g)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
    plot(yObserved(g).swe, yEstimated.(glacier), 'o', 'Color', options.RGB(g,:),'MarkerSize',2); hold on

        [F.(glacier), G.(glacier)] = fit(yObserved(g).swe, yEstimated.(glacier),'poly1');
        p = plot(F.(glacier)); hold on
            xlabel(''); ylabel('')
            set(p,'Color',options.RGB(g,:)); set(p, 'LineWidth',1.5);     
        xlabel('Observed WB (m w.e.)');
        if  g ==1;  ylabel('Estimated WB (m w.e.)');   end
                axis square;    box on;     grid on
        b = gca; legend(b,'off');
        title(options.GName{g})
        annotation('textbox',[locX(g) locY(g) .1 .1],...
                   'String', ['R^2=',num2str(round(G.(glacier).rsquare,2))],...
                   'EdgeColor','none','FontWeight','bold')
end

%saveFIG_IGS(['LRobservedVSestimated_',den],2,10)

%% Interp Method - Beta coeffs boxplot
    
    clear
load Full.mat fullLR 
load TestInterp.mat SigBetaDen
run OPTIONS.m

%Get betas for chosen density option
d = 2;
den = options.DenOpt{d};
    betas = fullLR.(den).coeff{1:7,:}; 
for g = 1:3;    glacier = options.glacier{g};
    betas(:,g) = betas(:,g).*SigBetaDen.(glacier)(2:end,d);
end

     incB = [1 2 4 6 7];
     
figure(8); clf        
B = bar(betas(incB,:));
        for g = 1:3;    glacier = options.glacier{g};
    B(g).FaceColor = options.RGB(g,:);
    B(g).EdgeColor = 'none';
        end
    legend(options.GName); % Add a legend
    ylabel('Regression coefficient'); hold on 
    F = gca; F.XTickLabel = {'Elevation','Distance from centreline','Slope',...
                             'Curvature','Wind redistribution'};
             F.XTickLabelRotation = 0;
        ylim([-0.05 0.11])
        for i = 1:4
        line([i+0.5 i+0.5],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end
 
saveFIG_IGS('BetaCoeffsBar',1,8)

%% WSMB Distribution - LR sources of var

%clear; load varWSMB.mat varB 
run OPTIONS.m

    close all
ylimmax = [50, 20, NaN, 50, 5];
legendX = [0.32, 0.60, 0, 0.32, 0.60];
legendY = [0.78, 0.78, 0, 0.31, 0.31];

figure(9); clf;
x = 0:0.001:1.2;
for o = 1:2
    if      o == 1; data = varB.LR.zz;      t = 'Grid scale';     
    elseif  o == 2; data = varB.LR.interp;  t = 'Interpolation';   
    end
for g = 1:3;     glacier = options.glacier{g};
for d = 1:8;     den = options.DenOpt{d};
    
ProbDen.(den).(glacier) = fitdist(data.(den).(glacier)(:),'Normal');
    y = pdf(ProbDen.(den).(glacier),x);  
    
subplot(1,3,o) 
p(g) = fill(x,y,options.RGB(g,:),'FaceAlpha',0.2, 'EdgeColor', 'none'); hold on
    if      o == 1;         ylabel({'Linear regression','probablity density'});    end
    xlabel('Glacier-wide WB (m w.e.)'); 
end
    ylim([0 ylimmax(o)]); 
    %legend
    if g == 3
        [L, icons] = legend(p,options.GName,'location','northeast'); 
        set(L, 'Position', [legendX(o), legendY(o), 0.01, 0.01])
            L.Box='off';
            for i = length(icons)/2+1:length(icons)
            icons(i).Vertices(3:4,1) = icons(i).Vertices(3:4,1)/2;
            end
            for i = 1:length(icons)/2
            icons(i).Position(1) = icons(i).Position(1)/2;
            end
    end

end
    xlim([min(x) 1])
    title(t);
end

%full PDF LR and SK

for g = 1:3
glacier = options.glacier{g};
for d = 1:8
den = options.DenOpt{d};
TLR.(glacier)(:,:,d)    = varB.LR.zzinterp.(den).(glacier);
end
end

for g = 1:3;     glacier = options.glacier{g};

ProbDenLR.(glacier) = fitdist(TLR.(glacier)(:),'Normal');
    yLR     = pdf(ProbDenLR.(glacier),x);   %yLR = yLR/Pmax(1);  
    std(yLR)
subplot(1,3,3) 
p(g) = fill(x,yLR,options.RGB(g,:),'FaceAlpha',0.8, 'EdgeColor', 'none'); hold on
    if g == 3
        [L1, icons1] = legend(p,options.GName,'location','northeast'); 
        set(L1, 'Position', [0.88, 0.78, 0.01, 0.01])
            L1.Box='off';
            for i = length(icons1)/2+1:length(icons1)
            icons1(i).Vertices(3:4,1) = icons1(i).Vertices(3:4,1)/2;
            end
            for i = 1:length(icons1)/2
            icons1(i).Position(1) = icons1(i).Position(1)/2;
            end
    end
        ylim([0 15])
        xlim([min(x) 1])
    title([{'Grid scale &'},{'Density assignment &'},{'Interpolation'}])
    xlabel('Glacier-wide WB (m w.e.)'); 

end


%saveFIG_IGS('WSMBDist',2,6)

%% WSMB Distribution - total spatial variability 

clear; load varWSMB.mat D* options
load TopoSWE.mat SWE

o = 1;    s = 'LR';   m = 0.6;
        
% Linear regression
for g = 1:3; glacier = options.glacier{g};
    DPlot.(glacier) = zeros(options.mapsize(g,:));
    for d = 1:8; den = options.DenOpt{d};
        DD.(glacier) = D.(s).(den).(glacier)/(max(D.(s).(den).(glacier)(:)));
        DD.(glacier)(options.mapNaN.(glacier)) = NaN;
    DPlot.(glacier) = DPlot.(glacier) + DD.(glacier);
    end
    DPlot.(glacier) = DPlot.(glacier)/nanmax(DPlot.(glacier)(:)*m);
end

figure(o);
PlotTopoParameter_IGS(DPlot,'summer','Relative uncertainity',SWE,'none','nomassB')
    saveFIG_IGS(['SpatialVar_',s],2,8.6)


%% Accumulation gradient
run OPTIONS
load Full.mat fullLR

taylor(:,1) = [571731.48;577258.98;580978.1;587346.4;591126.5;597353.2;601796.1;608101];
taylor(:,2) = [6737517.35;6733918.68;6730286.9;6730436.4;6724959.2;6730694.1;6734532.2;6736574.4];
taylor(:,3) = [2620;2640;2380;2225;2070;1915;1765;1615];
taylor(:,4) = [3.2;3.7;2.6;2.7;2.4;1.78;1.65;0.91];
taylor(:,5) = [0.407;0.409;0.394;0.364;0.388;0.385;0.390;0.342];
taylor(:,6) = [1.302;1.513;1.024;0.983;0.932;0.685;0.643;0.311];

alex(:,1) = [566453.4;570077.4;595349.9;601424.8;605031];
alex(:,2) = [6727621.2;6732429.7;6741065.2;6753607.6;6762717.2];
alex(:,3) = [2610;2730;2321;2472;2434];
alex(:,6) = [1.30;1.59;...
            nanmean(fullLR.S2.G4(:));...
            nanmean(fullLR.S2.G2(:));...
            nanmean(fullLR.S2.G13(:))];
    alexerr = [0.029,0.049,0.032];%*1.96;
    
    Dt = sqrt((taylor(:,1)-alex(1,1)).^2+(taylor(:,2)-alex(1,2)).^2)/1000;
    Da = sqrt((alex(:,1)-alex(1,1)).^2+(alex(:,2)-alex(1,2)).^2)/1000;
    [Dat, Iat] = sort([Da; Dt]); AT = [alex(:,6); taylor(:,6)]; AT = AT(Iat);    
[Fat, Gat] = fit(Dat,AT,'poly1');
    CIat = predint(Fat,Dat);    
    
figure(2); clf
    pt = plot(Fat);    set(pt,'Color','k'); set(pt, 'LineWidth',1.5); hold on
L(2) = plot(Da(1:2),alex(1:2,6),'o', 'MarkerSize',6, 'Color',[234, 140, 46]/255, 'MarkerFaceColor',[234, 140, 46]/255);
L(3) = errorbar(Da(3:5),alex(3:5,6),alexerr,'o', 'MarkerSize',6, 'Color',[234, 140, 46]/255);
L(1) = plot(Dt,taylor(:,6),'s', 'MarkerSize',6, 'Color',[68, 181, 226]/255, 'MarkerFaceColor',[68, 181, 226]/255); hold on
    xlim([min(Dat), max(Dat)+1]); 
    grid on
        xlabel('Distance from topographic divide (km)'); ylabel('WB (m w.e.)')
        LEG = legend(L,{'Taylor-Barge (1969)','Flowers unpublished (2016)','This study (2016)'});
            LEG.Position(1:2) = [0.58 0.755];
    ylim([0 1.7])
    
   %Labels 
    Lx = Da(3:5); Ly = alex(3:5,6);
    for g = 1:3
        strG = options.GName{g};
        if g == 1; bit = 0;
        else; bit = 0.1;
        end
        text(Lx(g,1)-1, Ly(g,1)+bit, strG,...
            'HorizontalAlignment','right','FontSize',9,'FontName','Arial');
    end
        
saveFIG_IGS('AccumGradPres',1,7.8)

