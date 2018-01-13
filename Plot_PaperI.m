%% Measurement - Depth Boxplot and SP vs FS
    clear
load TopoSWE.mat allDepth options Density

 %Depth boxplot
    lG4     = length(allDepth.G4(:)); 
    lG2     = length(allDepth.G2(:)); 
    lG13    = length(allDepth.G13(:));
    toBox = [[allDepth.G4(:); nan(lG13-lG4,1)],...
             [allDepth.G2(:); nan(lG13-lG2,1)],...
              allDepth.G13(:)];
figure(1); clf
subplot(2,1,1)
boxplot(toBox,'labels',{'Glacier 4','Glacier 2','Glacier 13'},...
               'BoxStyle','outline','Colors',options.RGB,'Width',0.35,...
               'OutlierSize',4,'Symbol','o')
    ylabel('Snow depth (cm)');
    annotation('textbox',[.15 .87 .1 .1],'String', 'a','EdgeColor','none','FontWeight','bold')


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

subplot(2,1,2)
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
    xlabel('SP-derived density (kg m^{-3})')
    ylabel('FS-derived density (kg m^{-3})')
     
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
        
    annotation('textbox',[.15 .39 .1 .1],'String', 'b','EdgeColor','none','FontWeight','bold')

saveFIG_IGS('DepthBoxplot_SPvsFS',1,15)

%% Grid Cell - Zigzag histrogram
    %clear; close all
load TopoSWE.mat SWEzz 
    
%     labels(:,1) = {'G4 LZ';'G4 MZ';'G4 UZ';'none'}; 
%     labels(:,2) = {'G2 LZ';'G2 MZ';'G2 UZ';'none'};
%     labels(:,3) = {'G13 LZ';'G13 M_1Z';'G13 M_2Z';'G13 UZ'};
    labels(:,1) = {'L';'M';'U';'none'}; 
    labels(:,2) = {'L';'M';'U';'none'};
    labels(:,3) = {'L';'M1';'M2';'U'};
    Fcolor      = [ 7,95,73;
                    13,191,149;
                    141,247,222;
                    162,135,2;
                    242,202,2;
                    254,234,134;
                    46,27,49
                    121, 74, 130;
                    214, 190, 218;
                    245, 239, 246]/255;
    STDzz = NaN(4,3);
figure(5); clf
c = 1;
    for g = 1:3
        zz = categories(SWEzz(g).ZZ);

    subplot(1,3,g)
        for j = 1:length(zz)
            clear ZZdata
            ZZdata = SWEzz(g).swe(SWEzz(g).ZZ==char(zz(j)));
            ZZdata = (ZZdata-mean(ZZdata));%/mean(ZZdata)*100;%/std(ZZdata);
% ZZdepth = SWEzz(g).depth(SWEzz(g).ZZ==char(zz(j)));
% STDzz(j,g) = std(ZZdepth);
                if j ==1;   bins    = 20;   %round(sqrt(length(SWEzz(g).swe)));
                            edges   = linspace(-0.15,0.15,bins); end
                            N       = histcounts(ZZdata,edges);   
            fill([min(edges) (edges(:,1:end-1)+edges(:,2:end))/2 max(edges)],[0 N/sum(N) 0],Fcolor(c,:),...
                  'EdgeColor','none','FaceAlpha',0.75); hold on 
            xlabel('b_w (m w.e.)');     ylabel('Probability')
            grid on
            xlim([-0.15 0.15])
            %title(options.glacier{g})
            %ax = gca; ax.XTick = [-40:20:40];
            c = c+1;
        end
            [L,icons] = legend(labels(1:length(zz),g));
                X = [0.30, 0.58, 0.85];
                Y = [0.63, 0.63, 0.615];
            L.Position(1:2) = [X(g), Y(g)];
            L.Position(3:4) = [0.05 0.3];
                L.Box='off';
                for i = length(icons)/2+1:length(icons)
                icons(i).Vertices(3:4,1) = icons(i).Vertices(3:4,1)/2;
                end
                for i = 1:length(icons)/2
                icons(i).Position(1) = icons(i).Position(1)/2;
                end
    end
            % ZZ maps
            Y = 0.6;    S = 0.32;
             %ZZmap = imread('/home/glaciology1/Documents/MastersDocuments/Paper I/ZZMapG4.jpeg');
             ZZmap = imread('/Users/Alexandra/Documents/SFU/MastersDocuments/Paper I/ZZMapG4.jpeg');
            axes('position',[0.062,Y+0.09,S*0.7,S*0.7]); 
            imshow(ZZmap);            axis off; 
             %ZZmap = imread('/home/glaciology1/Documents/MastersDocuments/Paper I/ZZMapG2.jpeg');
             ZZmap = imread('/Users/Alexandra/Documents/SFU/MastersDocuments/Paper I/ZZMapG2.jpeg');
            axes('position',[0.31,Y,S,S]); 
            imshow(ZZmap);            axis off;
             %ZZmap = imread('/home/glaciology1/Documents/MastersDocuments/Paper I/ZZMapG13.jpeg');
             ZZmap = imread('/Users/Alexandra/Documents/SFU/MastersDocuments/Paper I/ZZMapG13.jpeg');
            axes('position',[0.565,Y-0.025,S*1.1,S*1.1]); 
            imshow(ZZmap);            axis off;

    annotation('textbox',[.13 .95 .1 .1],'String', 'a','EdgeColor','none','FontWeight','bold')
    annotation('textbox',[.41 .95 .1 .1],'String', 'b','EdgeColor','none','FontWeight','bold')
    annotation('textbox',[.69 .95 .1 .1],'String', 'c','EdgeColor','none','FontWeight','bold')
            
            
saveFIG_IGS('ZigzagHistogram',2,6);

%% Interp Method - LR & SK map
    clear
load Full.mat fullLR fullOK options
load TopoSWE.mat SWE
for g = 1:3;    glacier = options.glacier{g};
    inputOK.(glacier) = fullOK.S2.(glacier).pred;
end

figure(6); clf
PlotTopoParameter_IGS(fullLR.S2, 'modelledSWE', 'b_w (m w.e.)', SWE, 'black', 'massB')
	saveFIG_IGS('LR_map',2,8.6)
figure(6); clf
PlotTopoParameter_IGS(inputOK, 'modelledSWE', 'b_w (m w.e.)', SWE, 'black', 'massB')
	saveFIG_IGS('OK_map',2,8.6)

%% Interp Method - Observed vs Estimated SWE
%     clear
% load Full.mat fullLR fullOK
% load TopoSWE.mat SWE topo_sampled options
OPTIONS

den = 'S2';
    yObserved   = ObsInCell(SWE, topo_sampled);
    for g = 1:3;    glacier = options.glacier{g};
    OKinput.(glacier) = fullOK.(den).(glacier).pred;
    end
    
          locX = [.06 .39 .72];     locX = [locX locX];
          locY = [.89; 0.41];       locY = repmat(locY,1,3); locY = [locY(1,:) locY(2,:)];

figure(7); clf
[ha, ~] = tight_subplot(2,3,[0.03 .01],[.05 0.01],[0.02 0]);
for y = 1:2
    if      y ==1;  yEstimated = SampledCell(fullLR.(den));   k = 0;
    elseif  y ==2;  yEstimated = SampledCell(OKinput);   k = 3;
    end

for g = 1:3;    glacier = options.glacier{g};
    
axes(ha(g+k))
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
    plot(yObserved(g).swe, yEstimated.(glacier), 'o', 'Color', options.RGB(g,:),'MarkerSize',2); hold on

        [F.(glacier), G.(glacier)] = fit(yObserved(g).swe, yEstimated.(glacier),'poly1');
        p = plot(F.(glacier)); hold on
            xlabel(''); ylabel('')
            set(p,'Color',options.RGB(g,:)); set(p, 'LineWidth',1.5);     
        if      y ==2; xlabel([{'Gridcell-averaged'}, {'WB (m w.e.)'}]); end
        if      y ==1 && g ==1;  ylabel([{'LR gridcell-'}, {'estimated b_w (m w.e.)'}]);
        elseif  y ==2 && g ==1;  ylabel([{'OK gridcell-'}, {'estimated b_w (m w.e.)'}]);
        end
                axis square;    box on;     grid on
        b = gca; legend(b,'off');
        
        annotation('textbox',[locX(g+k) locY(g+k) .1 .1],'String', options.glacier{g},'EdgeColor','none','FontWeight','bold')
        annotation('textbox',[locX(g+k) locY(g+k)-0.05 .1 .1],'String', ['R^2= ',num2str(round(G.(glacier).rsquare,2))],'EdgeColor','none')
end
end
            set(ha(1:6),'YTickLabel',num2cell(0:0.2:1)); set(ha(4:6),'XTickLabel',num2cell(0:0.5:1))


saveFIG_IGS(['observedVSestimated_',den],2,10)

%% Interp Method - Beta coeffs boxplot
    
    clear
load Full.mat fullLR 
load LR_SK_RK.mat SigBetaDen
OPTIONS

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

BETAS(:,:,[3,5]) = [];     
     incB = [1 2 4 6 7];
     
figure(8); clf        
aboxplot(BETAS,'labels',options.topoVars(incB), ...
         'Colormap', options.RGB, 'OutlierMarkerSize',8,...
         'WidthS',1.9,'WidthE',1.3); % Advanced box plot
        L = legend('G4','G2','G13'); % Add a legend
            set(L, 'Position', [0.7, 0.785, 0.05, 0.05])

            xlabel('Topographic parameters')
            ylabel('Regression coefficient'); hold on 
        
        ylim([-0.06 0.12])
        line(xlim,[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        for i = 1:4
        line([i+0.5 i+0.5],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end
 
saveFIG_IGS('BetaCoeffs',1,9)

%% WSMB Distribution - LR & SK sources of var

%clear; load varWSMB.mat varB options 

    close all
ylimmax = [50, 20, NaN, 50, 5];
legendX = [0.32, 0.60, 0, 0.32, 0.60];
legendY = [0.78, 0.78, 0, 0.31, 0.31];

figure(9); clf;
x = 0:0.001:1.2;
for o = [1,2,4,5]
    if      o == 1; data = varB.LR.zz;      t = '\sigma_{GS}';     
    elseif  o == 2; data = varB.LR.interp;  t = '\sigma_{INT}';   
    elseif  o == 4; data = varB.SK.zz;      t = '\sigma_{GS}';     
    elseif  o == 5; data = varB.SK.interp;  t = '\sigma_{INT}';    
    end
for g = 1:3;     glacier = options.glacier{g};
for d = 1:8;     den = options.DenOpt{d};
    
ProbDen.(den).(glacier) = fitdist(data.(den).(glacier)(:),'Normal');
    y = pdf(ProbDen.(den).(glacier),x);  
    
if  o == 5;    y = [0 y]; x = [x 0];    end

subplot(2,3,o) 
p(g) = fill(x,y,options.RGB(g,:),'FaceAlpha',0.2, 'EdgeColor', 'none'); hold on
    if      o == 1;         ylabel({'LR probablity','density'});
    elseif  o == 4;         ylabel({'SK probablity','density'});  
    end
    if  o == 4||o == 5; xlabel('Glacier-wide WB (m w.e.)');  end
end
    ylim([0 ylimmax(o)]); 
    %legend
    if g == 3
        [L, icons] = legend(p,options.glacier,'location','northeast'); 
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
TSK.(glacier)(:,:,d)    = varB.SK.zzinterp.(den).(glacier);
end
end

for g = 1:3;     glacier = options.glacier{g};

ProbDenLR.(glacier) = fitdist(TLR.(glacier)(:),'Normal');
    yLR     = pdf(ProbDenLR.(glacier),x);   %yLR = yLR/Pmax(1);  
ProbDenSK.(glacier) = fitdist(TSK.(glacier)(:),'Normal');
    ySK     = pdf(ProbDenSK.(glacier),x);   %ySK = ySK/Pmax(2);
    std(yLR)
subplot(2,3,3) 
p(g) = fill(x,yLR,options.RGB(g,:),'FaceAlpha',0.8, 'EdgeColor', 'none'); hold on
    if g == 3
        [L1, icons1] = legend(p,options.glacier,'location','northeast'); 
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
    title('\sigma_{GS} & \sigma_{\rho} & \sigma_{INT}')
subplot(2,3,6)     
q(g) = fill([x 0],[0 ySK],options.RGB(g,:),'FaceAlpha',0.8, 'EdgeColor', 'none'); hold on
    xlabel('B_w (m w.e.)'); 
    title('\sigma_{GS} & \sigma_{\rho} & \sigma_{INT}')
    xlim([min(x) 1])
    if g == 3
        [L, icons] = legend(q,options.glacier,'location','northeast'); 
        set(L, 'Position', [0.88, 0.31, 0.01, 0.01])
            L.Box='off';
            for i = length(icons)/2+1:length(icons)
            icons(i).Vertices(3:4,1) = icons(i).Vertices(3:4,1)/2;
            end
            for i = 1:length(icons)/2
            icons(i).Position(1) = icons(i).Position(1)/2;
            end
    end
end

saveFIG_IGS('WSMBDist',2,8)

%% WSMB Distribution - total spatial variability 

clear; load varWSMB.mat D* options
load TopoSWE.mat SWE

for o = 1:2
    if o == 1;      s = 'LR';   m = 0.6;
    elseif o == 2;  s = 'SK';   m = 0.95;
    end
        
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
end

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
%     plot(Dat, CIat(:,1),'--', 'Color','k', 'LineWidth',.2)
%     plot(Dat, CIat(:,2),'--','Color','k', 'LineWidth',.2)
L(2) = plot(Da(1:2),alex(1:2,6),'o', 'MarkerSize',6, 'Color',[234, 140, 46]/255, 'MarkerFaceColor',[234, 140, 46]/255);
L(3) = errorbar(Da(3:5),alex(3:5,6),alexerr,'o', 'MarkerSize',6, 'Color',[234, 140, 46]/255);
L(1) = plot(Dt,taylor(:,6),'s', 'MarkerSize',6, 'Color',[68, 181, 226]/255, 'MarkerFaceColor',[68, 181, 226]/255); hold on
    xlim([min(Dat), max(Dat)+1]); 
    ylim([0 2])
    grid on
        xlabel('Distance from topographic divide (km)'); ylabel('Winter balance (m w.e.)')
        LEG = legend(L,{'P-WB (1969)','P-WB (2016)','G-WB (2016)'});
            LEG.Position(1:2) = [0.58 0.755];
%             set(LEG,'PlotBoxAspectRatioMode','manual');
%             set(LEG,'PlotBoxAspectRatio',[1 0.8 1]);
    ylim([0 1.7])
    
   %Labels 
    Lx = Da(3:5); Ly = alex(3:5,6);
    for g = 1:3
        strG = options.glacier{g};
        if g == 1; bit = 0;
        else bit = 0.1;
        end
        text(Lx(g,1)-1, Ly(g,1)+bit, strG,...
            'HorizontalAlignment','right','FontSize',9,'FontName','Arial');
    end
        
saveFIG_IGS('AccumGrad',1,7.8)

%% SUPPLEMENTARY Topo Params - Full vs measured dist
clear
load TopoSWE.mat topo* options
    GT = {'Glacier 4','Glacier 2','Glacier 13'};

figure(1); clf
    header = fieldnames(topo_full.G4);
    n = 1;  nLab = 1:3:21;
for r = 1:3 %length(header)
    param = char(header(r));
    for i = 1:3
        name = char(options.glacier(i)); 
        %a(i) = subplot(length(header),3,n);
        a(i) = subplot(3,3,n);
            h(i).f = histogram(topo_full_ns.(name).(param)(:),25,'FaceColor','w'); hold on
            h(i).s = histogram(topo_sampled_ns.(name).(param),25,'FaceColor','k'); 
            
            hist_stats_full.(name)(r,:) = [nanmean(topo_full_ns.(name).(param)(:)),...
                                            nanstd(topo_full_ns.(name).(param)(:)),...
                                            skewness(topo_full_ns.(name).(param)(:)),...
                                            kurtosis(topo_full_ns.(name).(param)(:))];
            hist_stats_sample.(name)(r,:) = [nanmean(topo_sampled_ns.(name).(param)(:)),...
                                            nanstd(topo_sampled_ns.(name).(param)(:)),...
                                            skewness(topo_sampled_ns.(name).(param)(:)),...
                                            kurtosis(topo_sampled_ns.(name).(param)(:))];                            
                                        
            xlabel(char(options.topoVarsUnits(r)));     
            if ismember(n,nLab); ylabel('Frequency'); end
            axis tight
       %     legend('Full range','Sampled');
       if r ==1; title(GT(i)); end
    n = n+1;        
    end
    set([h(1).f h(1).s h(2).f h(2).s h(3).s], 'BinEdges', h(3).f.BinEdges)
    linkaxes(flip(a));
end 
saveFIG_IGS('TopoParamsSampled1',2,15)


figure(1); clf
n = 1;
for r = 4:7 %length(header)
    param = char(header(r));
    for i = 1:3
        name = char(options.glacier(i)); 
        a(i) = subplot(4,3,n);
            h(i).f = histogram(topo_full_ns.(name).(param)(:),25,'FaceColor','w'); hold on
            h(i).s = histogram(topo_sampled_ns.(name).(param),25,'FaceColor','k'); 
            
            hist_stats_full.(name)(r,:) = [nanmean(topo_full_ns.(name).(param)(:)),...
                                            nanstd(topo_full_ns.(name).(param)(:)),...
                                            skewness(topo_full_ns.(name).(param)(:)),...
                                            kurtosis(topo_full_ns.(name).(param)(:))];
            hist_stats_sample.(name)(r,:) = [nanmean(topo_sampled_ns.(name).(param)(:)),...
                                            nanstd(topo_sampled_ns.(name).(param)(:)),...
                                            skewness(topo_sampled_ns.(name).(param)(:)),...
                                            kurtosis(topo_sampled_ns.(name).(param)(:))];                            
                                        
            xlabel(char(options.topoVarsUnits(r)));     
            if ismember(n,nLab); ylabel('Frequency'); end
            axis tight
       %     legend('Full range','Sampled');
    n = n+1;        
    end
    set([h(1).f h(1).s h(2).f h(2).s h(3).s], 'BinEdges', h(3).f.BinEdges)
    linkaxes(flip(a));
end 
saveFIG_IGS('TopoParamsSampled2',2,20)