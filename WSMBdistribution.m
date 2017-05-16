%% WSMB distribution from fitlm coeffs
clear Q
varSIG = [0.027, 0.035, 0.04];

for g = 1:3
    glacier = options.glacier{g};
       
for d = 1:8
    den = options.DenOpt{d};

    varPD   = makedist('Normal','mu',0,'sigma',varSIG(g));
    varSWE  = random(varPD, length(fullSWE.(den).(glacier).swe),1000);
for vc = 1:1000
    dataSWE = fullSWE.(den).(glacier).swe + varSWE(:,vc);
    dataSWE(dataSWE<0) = 0;
    
 %Get Beta dist
    data    = [struct2table(topo_sampled.(glacier)), table(dataSWE,'VariableNames',{'swe'})];
    LRt     = fitlm(data);
    sigma   = LRt.CoefficientCovariance;
    mu      = LRt.Coefficients{:,1};
    beta    = mvnrnd(mu,sigma,1000);

for mc = 1:1000
     %Get random correlated beta value from normal distribution
    B = beta(mc,:);
    %B = LRt.Coefficients{:,1};
    
     %Get distribution of swe
    tempSWE.(glacier) = repmat(B(1),options.mapsize(g,:));

        fields = fieldnames(topo_full.(glacier));
    for f = 1:length(fields)
        param = fields{f};
        val = topo_full.(glacier).(param)*B(f+1);
        tempSWE.(glacier) = tempSWE.(glacier) + val; 
    end
        tempSWE.(glacier)(tempSWE.(glacier)<0) = 0;

 %Winter balance
Q.(den).(glacier)(mc,vc) = nanmean(tempSWE.(glacier)(:));
end
end
end
end
%% PLOT -> fitlm coeffs

 %all Gs, one density
bins    = 200;%round(sqrt(length(Qbeta.(den).G4(:))));
edges   = linspace(0.3,0.8,bins); 
for o = 1:3
    if      o == 1; data = Qzz;     t = '\sigma_{ZZ} Variability';      f = 'zz';
    elseif  o == 2; data = Qbeta;   t = '\beta Variability';             f = 'beta';
    elseif  o == 3; data = Qbetazz; t = '\beta and \sigma_{ZZ} Variability'; f = 'betaNzz';
    end
figure(4); clf; p = 1;
for d = 1:8
    den = options.DenOpt{d};
for g = 1:3
    glacier = options.glacier{g};
subplot(4,2,p)    
    h(g) = histogram(data.(den).(glacier)(:), edges, 'Normalization','probability',...
        'EdgeColor','none','FaceAlpha',0.5, 'FaceColor',options.RGB(g,:)); hold on
    ylabel('Probability'); xlabel('Winter surface mass balance (m w.e.)');  
end
    legend(h,options.glacier)
    title(den);
p = p+1;
end
saveFIG(['WSMB_Distribution',f],16)
end

    
 %With ZZ var vs Without
 c = [119, 24, 15; 71, 72, 137; 214, 164, 47]/255;
 figure(5); clf
 bins    = 200;%round(sqrt(length(Qbetazz.(den).G4(:))));
 edges   = linspace(0.3,0.8,bins); 
for g = 1:3
    glacier = options.glacier{g};
subplot(1,3,g)
    histogram(Qzz.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.6, 'FaceColor',c(1,:)); hold on
    histogram(Qbeta.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.6, 'FaceColor',c(2,:)); hold on
    histogram(Qbetazz.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.6, 'FaceColor',c(3,:)); hold on
   
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); title([glacier,' WSMB Distribution - ',den])
    legend('\sigma_{ZZ}','\beta','\beta and \sigma_{ZZ}','location','best') 
end
    saveFIG('WSMB_compareBetaAndZZBeta')
    
    
 %all Density, one G
figure(3); clf
c = cbrewer('qual','Set1',8);
for g = 1:3
glacier = options.glacier{g};

for d = 1:8
    den = options.DenOpt{d};

subplot(1,3,g)
    histogram(Qbeta.(den).(glacier), 'Normalization','probability',...
        'FaceColor',c(d,:),'EdgeColor','none','FaceAlpha',0.5); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); title(glacier)
    xlim([0.3 0.8])
end
    legend(options.DenOpt)
end
    saveFIG('WSMB_allDensityBeta')

 %all Density, one G
figure(2); clf
c = [147, 148, 150; 37, 37, 38]/255;   
p = 1;
for d = 1:8
    den = options.DenOpt{d};
for g = 1:3
glacier = options.glacier{g};

subplot(8,3,p)
    histogram(Qbeta.(den).(glacier), 'Normalization','probability',...
        'FaceColor',c(2,:),'EdgeColor','none'); hold on
    histogram(Qbetazz.(den).(glacier), 'Normalization','probability',...
        'FaceColor',c(1,:),'EdgeColor','none'); hold on
    ylabel('Prob.'); xlabel('WSMB (m w.e.)'); 
    if p <= 3; title([glacier, ' - ', den]); end
    p = p+1;
end
    legend('\beta','\beta and \sigma_{ZZ}')
end
    saveFIG('WSMB_allllll', 12)
    
    
%% Combining

for g = 1:3
glacier = options.glacier{g};
for d = 1:8
den = options.DenOpt{d};
Tbeta.(glacier)(:,:,d)  = Qbeta.(den).(glacier);
Tzz.(glacier)(:,:,d)    = Qzz.(den).(glacier);
Td.(glacier)(:,:,d)     = nanmean(fullLR.(den).(glacier));
Tdbetazz.(glacier)(:,:,d) = Qbetazz.(den).(glacier);
end
end


%dNbetaNzz
figure(1); clf
for g = 1:3
glacier = options.glacier{g};
histogram(Tdbetazz.(glacier)(:),'Normalization','probability',...
        'EdgeColor','none','FaceAlpha',0.75, 'FaceColor',options.RGB(g,:)); hold on
end
ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
legend(options.glacier)
        saveFIG('WSMB_dNbetaNzz')

%% WSMB from generous MLR coeffs

%load Full.mat fullCI
% for d = 2:9
%     den = options.DenOpt{d-1};
%     display(den)
% [ tempswe.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled); 
% 
% % Linear regression
% [~, fullCI.(den) ] =  LinearRegression( tempswe.(den), TOPOdata, topo_full );
%     
% end
clear Qfull
varSIG = [0.027, 0.035, 0.04];

for g = 1:3
    glacier = options.glacier{g};
       
for d = 1:8
    den = options.DenOpt{d};

%     varPD   = makedist('Normal','mu',0,'sigma',varSIG(g));
%     varSWE  = random(varPD, length(fullSWE.(den).(glacier).swe),1000);
for vc = 1%:1000
    dataSWE = fullSWE.(den).(glacier).swe;% + varSWE(:,vc);
    dataSWE(dataSWE<0) = 0;
    
 %Get Beta dist
    X       = [struct2table(topo_sampled.(glacier)), table(dataSWE,'VariableNames',{'swe'})]; X = X{:,:};
    sigma   = cov(X);
    mu      = fullLR.(den).coeff{[8,1:7],g};
    beta    = mvnrnd(mu,sigma,1000);

for mc = 1:1000
     %Get random correlated beta value from normal distribution
    B = beta(mc,:);
    %B = LRt.Coefficients{:,1};
    
     %Get distribution of swe
    tempSWE.(glacier) = repmat(B(1),options.mapsize(g,:));

        fields = fieldnames(topo_full.(glacier));
    for f = 1:length(fields)
        param = fields{f};
        val = topo_full.(glacier).(param)*B(f+1);
        tempSWE.(glacier) = tempSWE.(glacier) + val; 
    end
        tempSWE.(glacier)(tempSWE.(glacier)<0) = 0;

 %Winter balance
fullQ.(den).(glacier)(mc,vc) = nanmean(tempSWE.(glacier)(:));
end
end
end
end

%% PLOT -> full generous coeffs

 %all Gs, one density
den = options.DenOpt{7};

figure(4); clf
for g = 1:3
    glacier = options.glacier{g};

%subplot(1,3,g)
    histogram(fullQ.(den).(glacier), 'EdgeColor','none','FaceAlpha',0.7); hold on
    ylabel('Frequency'); xlabel('Winter surface mass balance'); title(glacier)
end
    legend(options.glacier)
    
 %all Density, one G
for g = 1:3;
glacier = options.glacier{g};

figure(g); clf
for d = 1:8
    den = options.DenOpt{d};

%subplot(1,3,g)
    histogram(fullQ.(den).(glacier), 'EdgeColor','none','FaceAlpha',0.7); hold on
    ylabel('Frequency'); xlabel('Winter surface mass balance'); title(glacier)
end
    legend(options.DenOpt)
end


%% KRIGING

clear Q
load TopoSWE.mat
varSIG = [0.027, 0.035, 0.04];
    SWE = ObsInCell(SWE, topo_sampled);
 for g = 1:3
        glacier = options.glacier{g};
            varPD.(glacier)   = makedist('Normal','mu',0,'sigma',varSIG(g));
            varSWE.(glacier)  = random(varPD.(glacier), length(SWE(g).swe),1000);
 end


for d = 1:8
    den = options.DenOpt{d};
    display(den)
for vc = 1:1000
        [ dataSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled); 
    for gg = 1:3
        glacier = options.glacier{gg};           
            dataSWE.(den).(glacier)(:,1) = dataSWE.(den).(glacier)(:,1) + varSWE.(glacier)(:,vc);
                dataSWE.(den).(glacier)(dataSWE.(den).(glacier)(:,1)<0,1) = 0;
    end
    
% Simple kriging

  tempSK.(den) =  KrigingR_G( dataSWE.(den) );

 %Winter balance
for gg = 1:3
        glacier = options.glacier{gg};      
        Q.(den).(glacier)(vc) = nanmean(tempSK.(den).(glacier)(:));
end
end
end
