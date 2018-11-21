%% bw Map

%load Full.mat fullLR fullOK options
load TopoSWE.mat SWE
for g = 1:3;    glacier = options.glacier{g};
    inputOK.(glacier) = fullOK.S2.(glacier).pred;
end

figure(6); clf
PlotTopoParameter_IGS(inputOK, 'modelledSWE', 'b_w (m w.e.)', SWE, 'black', 'massB')
% 	saveFIG_IGS('OK_map',2,8.6)

%% STD MAP

 % INTERPOLATION
%OK     
load Full.mat
for g = 1:3;    glacier = options.glacier{g};
    stdMap.OK.int.(glacier) = fullOK.S2.(glacier).std;
    max(stdMap.OK.int.(glacier)(:))
end

    
 % DENSITY
 for i = 1:2
    for g = 1:3;        glacier = options.glacier{g};
        stack = zeros(options.mapsize(g,1),options.mapsize(g,2),8);
    for d = 1:8;        den = options.DenOpt{d};
        if i == 1;      stack(:,:,d) = fullLR.(den).(glacier);        t = 'LR'; 
        elseif i ==2;   stack(:,:,d) = fullSK.(den).(glacier).pred;   t = 'OK';
        end
    end
        stdMap.(t).den.(glacier) = std(stack,[],3);
            max(stdMap.(t).den.(glacier)(:))
    end
  end
 
 % GRIDCELL
 %LR
        %load('LRstdGS.mat')
     for g = 1:3;        glacier = options.glacier{g};
         stack = zeros(options.mapsize(g,1),options.mapsize(g,2),1000);
     for mc = 1:1000
        stack(:,:,mc) = LRstdGS(mc).(glacier);
     end
     stdMap.LR.gs.(glacier) = std(stack,[],3);
     end
     
 %OK    
     for g = 1:3;        glacier = options.glacier{g};
         stack = zeros(options.mapsize(g,1),options.mapsize(g,2),500);
     for mc = 1:500
        stack(:,:,mc) = KRIGzz.S2(mc).(glacier).pred;
     end
     stdMap.OK.gs.(glacier) = std(stack,[],3);
     end


% COMBINING ALL
     for g = 1:3;        glacier = options.glacier{g};
         clear S
        S(:,:,1) = stdMap.LR.int.(glacier);
        S(:,:,2) = stdMap.LR.den.(glacier); 
        S(:,:,3) = stdMap.LR.gs.(glacier);
        stdMap.LR.ALL.(glacier) = sqrt(sum(S.^2,3));
        
        S(:,:,1) = stdMap.OK.int.(glacier);
        S(:,:,2) = stdMap.OK.den.(glacier); 
        S(:,:,3) = stdMap.OK.gs.(glacier);
        stdMap.OK.ALL.(glacier) = sqrt(sum(S.^2,3));
            max(stdMap.LR.ALL.(glacier)(:))
     end
   %%  
   load stdMap.mat
   lr_scale = 0.1;
   ok_scale = 0.2;
   
figure(6); clf
PlotTopoParameter_IGS(stdMap.LR.gs, 'summer', 'STD of b_w (m w.e.)', SWE, 'black', 'massB',lr_scale)
	saveFIG_IGS('LRstd_map_gridcell',2,8.6)
figure(7); clf
PlotTopoParameter_IGS(stdMap.OK.gs, 'summer', 'STD of b_w (m w.e.)', SWE, 'black', 'massB',ok_scale)
	saveFIG_IGS('OKstd_map_gridcell',2,8.6)
     
figure(6); clf
PlotTopoParameter_IGS(stdMap.LR.den, 'summer', 'STD of b_w (m w.e.)', SWE, 'black', 'massB',lr_scale)
	saveFIG_IGS('LRstd_map_density',2,8.6)
figure(7); clf
PlotTopoParameter_IGS(stdMap.OK.den, 'summer', 'STD of b_w (m w.e.)', SWE, 'black', 'massB',ok_scale)
	saveFIG_IGS('OKstd_map_density',2,8.6)
  
figure(6); clf
PlotTopoParameter_IGS(stdMap.LR.int, 'summer', 'STD of b_w (m w.e.)', SWE, 'black', 'massB',lr_scale)
	saveFIG_IGS('LRstd_map_int',2,8.6)
figure(7); clf
PlotTopoParameter_IGS(stdMap.OK.int, 'summer', 'STD of b_w (m w.e.)', SWE, 'black', 'massB',ok_scale)
	saveFIG_IGS('OKstd_map_int',2,8.6)
         
figure(6); clf
PlotTopoParameter_IGS(stdMap.LR.ALL, 'summer', 'Combined STD of b_w (m w.e.)', SWE, 'black', 'massB',lr_scale)
	saveFIG_IGS('LRstd_map_ALL',2,8.6)
figure(7); clf
PlotTopoParameter_IGS(stdMap.OK.ALL, 'summer', 'Combined STD of b_w (m w.e.)', SWE, 'black', 'massB',ok_scale)
	saveFIG_IGS('OKstd_map_ALL',2,8.6)

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

    ylimmax = [50, 20, 15, 240, 5, 5]; %Manual control over y limits on plot (feel free to change them)

    % Legend locations    
legendX = [0.27, 0.55, 0.83, 0.27, 0.55, 0.83];
legendY = [0.74, 0.74, 0.74, 0.26, 0.26, 0.26];

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
    elseif  o == 4;         ylabel({'OK probablity','density'});  
    end
    if  o == 4||o == 5; xlabel('Glacier-wide WB (m w.e.)');  end
end
    ylim([0 ylimmax(o)]); 
    %legend
    if g == 3
        [L, icons] = legend(p,options.glacier,'location','northeast'); 
        set(L, 'Position', [legendX(o), legendY(o), 0.09, 0.16])
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
        set(L, 'Position', [legendX(o), legendY(o), 0.09, 0.16])
            L.Box='off';
            for i = length(icons)/2+1:length(icons)
            icons(i).Vertices(3:4,1) = icons(i).Vertices(3:4,1)/2;
            end
            for i = 1:length(icons)/2
            icons(i).Position(1) = icons(i).Position(1)/2;
            end
    end

    ylim([0 ylimmax(o)])
    xlim([min(x) 1])
    title('\sigma_{GS} & \sigma_{\rho} & \sigma_{INT}')
    if  o == 6; xlabel('Glacier-wide WB (m w.e.)');  end

end
end
    
saveFIG_IGS('WSMBDist',2,9.5)

%% WSMB Distribution - total spatial variability 

% load varWSMB.mat D* 
% load MonteCarloOK.mat DOK
 load TopoSWE.mat SWE
% run OPTIONS
    D.OK = DOK;

for o = 2%1:2
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
    saveFIG_IGS(['SpatialVar_',s],2,8.6)
end


%% 1000 runs vs 500 runs (S1)

for g = 1:3;    glacier = options.glacier{g};
for d = 1;    den = options.DenOpt{d};

%Calculate glacier-wide winter balance and std for all runs
%(***sigms_GS***)
    for mc = 1:500
        BwKRIGzz500.(den).(glacier)(mc) = nanmean(KRIGzz500.(den)(mc).(glacier).pred(:));               
    end
    for mc = 1:1000
        BwKRIGzz1000.(den).(glacier)(mc) = nanmean(KRIGzz1000.(den)(mc).(glacier).pred(:));               
    end
end
end
        clf;
        M = zeros(2,3); S = M;
for g = 1:3; glacier = options.glacier{g};
   subplot(1,3,g)
   histogram(BwKRIGzz1000.S1.(glacier),20,'FaceAlpha',0.7,'EdgeColor','none'); hold on
   histogram(BwKRIGzz500.S1.(glacier),20,'FaceAlpha',0.7,'EdgeColor','none')
        if g ==3; legend('1000 runs','500 runs'); end
        title(glacier)
        xlabel('Glacier-wide WB (m w.e.)'); 
        if g ==1; ylabel({'OK probablity','density'}); end
   M(1,g) = mean(BwKRIGzz500.S1.(glacier));        M(2,g) = mean(BwKRIGzz1000.S1.(glacier));
   S(1,g) = std(BwKRIGzz500.S1.(glacier));         S(2,g) = std(BwKRIGzz1000.S1.(glacier));

end
    saveFIG_IGS('MCruns1000vs500',2,8.6)
%% Updating tables

%Table 3 (Bw and RMSE)
%OKsigmaALL{1,:}

load TopoSWE.mat 
load Full.mat
    yObserved  = ObsInCell(SWE, topo_sampled);

        RMSE = zeros(8,3);
for d = 1:8;    den = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
%     InP.(glacier) = fullOK.(den).(glacier).pred;
    InP.(glacier) = fullLR.(den).(glacier);
end
    est.(den) = SampledCell(InP); %Model values
end

for d = 1:8;    den = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
    RMSE(d,g) = sqrt(mean((est.(den).(glacier)-yObserved(g).swe).^2));
end
end
mean(RMSE)


for d = 3;    den = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
%     Bw = nanmean(nanmean(fullOK.(den).(glacier).pred))
    Bw = nanmean(nanmean(fullLR.(den).(glacier)))
end
end


% % mean(RMSE)./OKsigmaALL{1,:}*100
% % 
% % 
% % 
% % %Table 4 (std of distribution)
% % mean(OKsigmaGS.std{:,:})*10^2
% % OKsigmaRHO{2,:}*10^2
% % mean(OKsigmaINT.std{:,:})*10^2
% % OKsigmaALL{2,:}*10^2