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

fig = figure(1);
for i = 1:length(SP)
errorbarxy(SP(i),FS(i),errorSP(i,2),errorFS(i,2),errorSP(i,1),errorFS(i,1),...
                'Color','k','LineStyle','none','Marker',markerS{i},...
                'MarkerFaceColor',markerC(i,:),'LineWidth',1,'MarkerSize',15,...
                'MarkerEdgeColor','none'); hold on
end
    axis([220 400 220 400])
    grid on
    line = refline(1,0);
        line.Color = 'k'; line.LineStyle = '--'; hold on
    xlabel('Snow-pit-derived density (kg m^{-3})')
    ylabel({'Federal-Sampler-','derived density (kg m^{-3})'})
     
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
%     fig.PaperUnits = 'centimeters';
%     fig.PaperSize  = [1000 1000];

        
saveFIG('SPvsFS',24)


%% Interp Method - LR map
    clear
load Full.mat fullLR 
load TopoSWE.mat SWE

PlotTopoParameter(fullLR.S2, 'modelledSWE', 'WB (m w.e.)', SWE, 'black', 'massB')
	saveFIG('LR_map',20)


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
    
          locX = [.14 .43 .71];     
          locY = 0.73;       

figure(7); clf
yEstimated = SampledCell(fullLR.(den)); 
   
for g = 1:3;    glacier = options.glacier{g};
    
subplot(1,3,g)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
    plot(yObserved(g).swe, yEstimated.(glacier), 'o', 'Color', options.RGB(g,:),'MarkerSize',5); hold on

        [F.(glacier), G.(glacier)] = fit(yObserved(g).swe, yEstimated.(glacier),'poly1');
        p = plot(F.(glacier)); hold on
            xlabel(''); ylabel('')
            set(p,'Color',options.RGB(g,:)); set(p, 'LineWidth',2.5);     
        xlabel('Observed WB (m w.e.)');
        if  g ==1;  ylabel('Estimated WB (m w.e.)');   end
                axis square;    box on;     grid on
        b = gca; legend(b,'off');
        title(options.GName{g})
        annotation('textbox',[locX(g) locY .1 .1],...
                   'String', ['R^2=',num2str(round(G.(glacier).rsquare,2))],...
                   'EdgeColor','none','FontWeight','bold')
end

saveFIG(['LRobservedVSestimated_',den],20)

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
    F = gca; F.XTickLabel = {'Elevation','Centre dist.','Slope',...
                             'Curvature','Wind redist.'};
             F.XTickLabelRotation = 50;
        ylim([-0.05 0.11])
        for i = 1:4
        line([i+0.5 i+0.5],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end
 
saveFIG('BetaCoeffsBar',20)

%% WSMB Distribution - LR sources of var

%clear; load varWSMB.mat varB 
run OPTIONS.m

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
for d = 1:8;     
    den = options.DenOpt{d};
    
ProbDen.(den).(glacier) = fitdist(data.(den).(glacier)(:),'Normal');
    y = pdf(ProbDen.(den).(glacier),x);  
    
subplot(1,3,o) 
p(g) = fill(x,y,options.RGB(g,:),'FaceAlpha',0.4, 'EdgeColor', 'none'); hold on
    if      o == 1;         ylabel({'Linear regression','probablity density'});    end
    xlabel('Glacier-wide WB (m w.e.)'); 
end
    ylim([0 ylimmax(o)]); 
    %legend
%     if g == 3
%         [L, icons] = legend(p,options.GName,'location','northeast'); 
%         set(L, 'Position', [legendX(o), legendY(o), 0.01, 0.01])
%             L.Box='off';
%             for i = length(icons)/2+1:length(icons)
%             icons(i).Vertices(3:4,1) = icons(i).Vertices(3:4,1)/2;
%             end
%             for i = 1:length(icons)/2
%             icons(i).Position(1) = icons(i).Position(1)/2;
%             end
%     end

end
    xlim([0.25 0.75])

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
%     if g == 3
%         [L1, icons1] = legend(p,options.GName,'location','northeast'); 
%         set(L1, 'Position', [0.88, 0.78, 0.01, 0.01])
%             L1.Box='off';
%             for i = length(icons1)/2+1:length(icons1)
%             icons1(i).Vertices(3:4,1) = icons1(i).Vertices(3:4,1)/2;
%             end
%             for i = 1:length(icons1)/2
%             icons1(i).Position(1) = icons1(i).Position(1)/2;
%             end
%     end
        ylim([0 15])
        xlim([0.25 0.75])
    title([{'Depth & Density & Interpolation'}])
    xlabel('Glacier-wide WB (m w.e.)'); 

end


saveFIG('WSMBDistLR',16)

% figure(2); clf
% for g = 1:3;     glacier = options.glacier{g};
% ProbDenLR.(glacier) = fitdist(TLR.(glacier)(:),'Normal');
%     yLR     = pdf(ProbDenLR.(glacier),x);   %yLR = yLR/Pmax(1);  
%     std(yLR)
% subplot(1,3,g) 
% p(g) = fill(x,yLR,options.RGB(g,:),'FaceAlpha',0.8, 'EdgeColor', 'none'); hold on
%         ylim([0 15])
%         xlim([0.25 0.75])
%         xlabel('');         ylabel(''); 
%         %set(gca,'xticklabel',{[]}) 
%         set(gca,'yticklabel',{[]}) 
% end
% saveFIG('WSMBDisttemp',16)

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
    
    Msize = 16;
    
figure(2); clf
    pt = plot(Fat);    set(pt,'Color','k'); set(pt, 'LineWidth',2); hold on
L(2) = plot(Da(1:2),alex(1:2,6),'^', 'MarkerSize',Msize, 'Color',[234, 140, 46]/255, 'MarkerFaceColor',[234, 140, 46]/255);
L(3) = errorbar(Da(3:5),alex(3:5,6),alexerr,'o', 'MarkerSize',Msize-6, 'Color',[24, 145, 26]/255, 'MarkerFaceColor',[24, 145, 26]/255);
L(3).LineWidth = 2;
L(1) = plot(Dt,taylor(:,6),'s', 'MarkerSize',Msize, 'Color',[68, 181, 226]/255, 'MarkerFaceColor',[68, 181, 226]/255); hold on
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
        if g == 1; bit = 0; bx = 1;
        else bit = 0.1; bx = -1;
        end
        text(Lx(g,1)-bx, Ly(g,1)+bit, strG,...
            'HorizontalAlignment','right','FontName','Arial');
    end
        
saveFIG('AccumGradPres',22)

%% G13 basic sampling

load Patterns.mat
clf
lw = 6;
Eg = options.rig.G13(:,1); Ng = options.rig.G13(:,2); 
plot(Eg(1:304),Ng(1:304),'k','LineWidth',lw); hold on
plot(Eg(305:330),Ng(305:330),'k','LineWidth',lw);
plot(Eg(331:356),Ng(331:356),'k','LineWidth',lw);
plot(Eg(357:end),Ng(357:end),'k','LineWidth',lw);
plot(pUTM.Centreline.G13(1:20:120,1),pUTM.Centreline.G13(1:20:120,2),'.','MarkerSize',70,'Color',[45, 150, 48]/255) %summer
%plot(pUTM.Centreline.G13(1:20:120,1),pUTM.Centreline.G13(1:20:120,2),'.','MarkerSize',70,'Color',[69 129 142]/255) %winter
axis off
axis square

saveFIG('G13TypicalS',22)



clf; load('Full.mat')
data = fullLR.S2.G13;
h = imagesc(data); hold on
set(h,'alphadata',~isnan(data));
 C = cbrewer('div','RdBu',40);
 %colormap(C);
 colormap(C(25:40,:));
%c = colorbar;   
axis square
axis off
saveFIG('G13TypicalWB',22)
