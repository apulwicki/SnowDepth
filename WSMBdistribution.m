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
den = options.DenOpt{7};

bins    = 200;%round(sqrt(length(Qbeta.(den).G4(:))));
edges   = linspace(0.3,0.8,bins); 
for d = 1:3
    if      d == 1; data = Qzz;     t = '\sigma_{ZZ} Variability';      f = 'zz';
    elseif  d == 2; data = Qbeta;   t = '\beta Variability';             f = 'beta';
    elseif  d == 3; data = Qbetazz; t = '\beta and \sigma_{ZZ} Variability'; f = 'betaNzz';
    end
figure(4); clf
for g = 1:3
    glacier = options.glacier{g};
    
    histogram(data.(den).(glacier), edges, 'Normalization','probability',...
        'EdgeColor','none','FaceAlpha',0.7, 'FaceColor',options.RGB(g,:)); hold on
    ylabel('Probability'); xlabel('Winter surface mass balance (m w.e.)'); 
    title({['WSMB Distribution - ',den],t})
end
    legend(options.glacier)
    saveFIG(['WSMB_',f])
end
    
 %With ZZ var vs Without
 c = [0, 53, 140; 57, 103, 178; 145, 175, 224]/255;
 figure(5); clf
 bins    = round(sqrt(length(Qbetazz.(den).G4(:))));
 edges   = linspace(0,1,bins); 
for g = 1:3
    glacier = options.glacier{g};
subplot(1,3,g)
    histogram(Qzz.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.7, 'FaceColor',c(1,:)); hold on
    histogram(Qbeta.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.7, 'FaceColor',c(2,:)); hold on
    histogram(Qbetazz.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.7, 'FaceColor',c(1,:)); hold on
   
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); title([glacier,' WSMB Distribution - ',den])
    legend('\sigma_{ZZ}','\beta','\beta and \sigma_{ZZ}') 
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
    histogram(Qbetazz.(den).(glacier), 'Normalization','probability',...
        'FaceColor',c(d,:),'EdgeColor','none','FaceAlpha',0.5); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); title(glacier)
    xlim([0 1.7])
    %xlim([0.3 0.8])
end
    legend(options.DenOpt)
end
    saveFIG('WSMB_allDensityBetaNzz')

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
%% WSMB from generous MLR coeffs

load Full.mat fullCI
% for d = 2:9
%     den = options.DenOpt{d-1};
%     display(den)
% [ tempswe.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled); 
% 
% % Linear regression
% [~, fullCI.(den) ] =  LinearRegression( tempswe.(den), TOPOdata, topo_full );
%     
% end

for g = 1:3
    glacier = options.glacier{g};

for d = 1:8
    den = options.DenOpt{d};
    
 %Get CI
SGt = diff(fullCI.(den).(glacier),1,2)/5.56;
MEt = mean(fullCI.(den).(glacier),2);
for i = 1:length(fullCI.(den).(glacier));
    pd(i) = makedist('Normal','mu',MEt(i),'sigma',SGt(i));
end

for mc = 1:10000;
 %Get random beta value from normal distribution
    B = zeros(size(MEt));
for i = 1:length(fullCI.(den).(glacier));
    B(i) = random(pd(i));
end
% figure(1); plot(pd); hold on; legend(LRt.CoefficientNames);

 %Get distribution of swe
tempSWE.(glacier) = repmat(B(end),options.mapsize(g,:));

    fields = fieldnames(topo_full.(glacier));
for f = 1:length(fields)
    param = fields{f};
    val = topo_full.(glacier).(param)*B(f);
    tempSWE.(glacier) = tempSWE.(glacier) + val; 
end
    tempSWE.(glacier)(tempSWE.(glacier)<0) = 0;
    
 %Winter balance
fullQ.(den).(glacier)(mc,1) = nanmean(tempSWE.(glacier)(:));

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