%% WSMB Distribution - LR & SK sources of var

%clear; close all
load varWSMB.mat varB options %varB is large so this might take a while
load MonteCarloOK.mat
    

    ylimmax = [50, 20, NaN, 50, 5]; %Manual control over y limits on plot (feel free to change them)

    % Legend locations    
legendX = [0.32, 0.60, 0, 0.32, 0.60];
legendY = [0.78, 0.78, 0, 0.31, 0.31];

figure(9); clf;
x = 0:0.001:1.2;
for o = [1,2]
    if      o == 1; data = varB.LR.zz;      t = '\sigma_{GS}';     
    elseif  o == 2; data = varB.LR.interp;  t = '\sigma_{INT}';   
    elseif  o == 4; data = BwKRIGzz;        t = '\sigma_{GS}';     
    elseif  o == 5; data = BwKRIGinterp;    t = '\sigma_{INT}';    
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
TSK.(glacier)           = BwKRIGall.(glacier);
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

load varWSMB.mat D* 
load MonteCarloOK.mat DOK
load TopoSWE.mat SWE
run OPTIONS

for o = 1:2
    if o == 1;      s = 'LR';   m = 0.6; %Change these m values to scale the color accordinly
    elseif o == 2;  s = 'OK';   m = 0.95;
    end
        
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
