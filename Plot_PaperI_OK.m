%% WSMB Distribution - LR & SK sources of var

%clear; close all
%load varWSMB.mat varB options %varB is large so this might take a while
%load MonteCarloOK.mat
    
%Fitting normal distributions to LR results (obtain mean and std for
%plotting)
        LRsigmaGS.mean  = table(zeros(8,1),zeros(8,1),zeros(8,1),'VariableNames',options.glacier,'RowNames',options.DenOpt);
        LRsigmaGS.std   = LRsigmaGS.mean;
        LRsigmaINT.mean = LRsigmaGS.mean;    LRsigmaINT.std = LRsigmaGS.mean;
for g = 1:3;     glacier = options.glacier{g};
for d = 1:8;     den = options.DenOpt{d};

T = fitdist(varB.LR.zz.(den).(glacier)(:),'Normal');
        LRsigmaGS.mean{d,g} = T.mu;     LRsigmaGS.std{d,g} = T.sigma;             
T = fitdist(varB.LR.interp.(den).(glacier)(:),'Normal');
        LRsigmaINT.mean{d,g} = T.mu;    LRsigmaINT.std{d,g} = T.sigma;             

end
end

    ylimmax = [50, 20, 15, 180, 5, 5]; %Manual control over y limits on plot (feel free to change them)

    % Legend locations    
legendX = [0.32, 0.60, 0.90, 0.32, 0.60, 0.90];
legendY = [0.78, 0.78, 0.78, 0.31, 0.31, 0.31];

x = 0:0.001:1.2;

figure(9); clf;
for o = [1,2,4,5]
    Falpha = 0.2;
    if      o == 1;     t = '\sigma_{GS}';  
        mu = LRsigmaGS.mean;    sigma = LRsigmaGS.std;        
    elseif  o == 2;     t = '\sigma_{INT}';
        mu = LRsigmaINT.mean;   sigma = LRsigmaINT.std;        
    elseif  o == 4;     t = '\sigma_{GS}';    Falpha = 0.7;  
        mu = OKsigmaGS.mean;    sigma = OKsigmaGS.std;        
    elseif  o == 5;     t = '\sigma_{INT}'; 
        mu = OKsigmaINT.mean;   sigma = OKsigmaINT.std;        
    end
for g = 1:3;     glacier = options.glacier{g};
for d = 1:8;     den = options.DenOpt{d};
    
    y = normpdf(x, mu{d,g}, sigma{d,g});  
    
subplot(2,3,o); 
p(g) = fill([x 0],[0 y],options.RGB(g,:),'FaceAlpha',Falpha, 'EdgeColor', 'none'); hold on
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

for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};
    TLR.(glacier)(:,:,d)    = varB.LR.zzinterp.(den).(glacier);
    TSK.(glacier)           = BwKRIGall.(glacier);
end
end

for o = [3,6]
for g = 1:3;     glacier = options.glacier{g};
    if      o == 3;     t = '\sigma_{GS}';  
        ProbDenLR.(glacier) = fitdist(TLR.(glacier)(:),'Normal');
    y = pdf(ProbDenLR.(glacier),x);   %yLR = yLR/Pmax(1);  
    elseif  o == 6;     t = '\sigma_{INT}';
    y = normpdf(x, OKsigmaALL{1,g}, OKsigmaALL{2,g});   %ySK = ySK/Pmax(2);
    end

subplot(2,3,o) 
p(g) = fill([x 0],[0 y],options.RGB(g,:),'FaceAlpha',0.8, 'EdgeColor', 'none'); hold on
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
%     if g == 3
%         [L1, icons1] = legend(p,options.glacier,'location','northeast'); 
%         set(L1, 'Position', [0.88, 0.78, 0.01, 0.01])
%             L1.Box='off';
%             for i = length(icons1)/2+1:length(icons1)
%             icons1(i).Vertices(3:4,1) = icons1(i).Vertices(3:4,1)/2;
%             end
%             for i = 1:length(icons1)/2
%             icons1(i).Position(1) = icons1(i).Position(1)/2;
%             end
%     end
    ylim([0 ylimmax(o)])
    xlim([min(x) 1])
    title('\sigma_{GS} & \sigma_{\rho} & \sigma_{INT}')
    if  o == 6; xlabel('Glacier-wide WB (m w.e.)');  end

end
end
    
saveFIG_IGS('WSMBDist',2,8)

%% WSMB Distribution - total spatial variability 

load varWSMB.mat D* 
load MonteCarloOK.mat DOK
load TopoSWE.mat SWE
run OPTIONS
    D.OK = DOK;

for o = 1:2
    if o == 1;      s = 'LR';   m = 0.6; %Change these m values to scale the color accordinly
    elseif o == 2;  s = 'OK';   m = 0.6;
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
%    saveFIG_IGS(['SpatialVar_',s],2,8.6)
end
