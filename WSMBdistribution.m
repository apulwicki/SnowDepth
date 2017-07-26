
%% WSMB - SWE Var

format shortg
clock
t = cputime;
%load Full.mat; load TopoSWE.mat
%clear fullQzz

for d = 1:8
    den = options.DenOpt{d};
[ inputSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled);
end

varSIG = [0.027, 0.035, 0.04];
for g = 1:3;    glacier = options.glacier{g};
    varPD.(glacier)   = makedist('Normal','mu',0,'sigma',varSIG(g));
    varSWE.(glacier)  = random(varPD.(glacier), length(inputSWE.S1.(glacier)),1000);   
end

mc = 1;

for d = 1:4
    den = options.DenOpt{d};
    display(den)
for vc = 1:1000
    for g = 1:3;        glacier = options.glacier{g};
        dataSWE.(den).(glacier)(:,1) = inputSWE.(den).(glacier)(:,1) + varSWE.(glacier)(:,vc);
        I = (dataSWE.(den).(glacier)(:,1)<0);
        dataSWE.(den).(glacier)(I,1) = 0;
    end
   
% Linear regression

    LRtemp.(den)(vc)  =  LinearRegression( dataSWE.(den), TOPOdata, topo_full );
    
    for g = 1:3;        glacier = options.glacier{g};
 %Winter balance
fullQzz.(den).(glacier)(vc) = nanmean(LRtemp.(den)(vc).(glacier)(:));
    end
    
end
clock
t = (cputime-t)/60/60
end



%% WSMB - Beta Var

format shortg
clock
t = cputime;
%load Full.mat; load TopoSWE.mat
%clear fullQbeta

for d = 1:8
    den = options.DenOpt{d};
[ dataSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled);
    
    for g = 1:3;        glacier = options.glacier{g};
        mu      = fullLR.(den).coeff{[8,1:7],g};
    
 %Get Beta dist
    X       = struct2table(TOPOdata.(glacier)); X = X{:,:}; X = [ones(length(X),1) X];
    Y       = dataSWE.(den).(glacier)(:,1);
    
    %VarCov matrix from betas
    sigmaSq = sum((Y-X*mu).^2)/(size(X,1)-size(X,2));
    VarCov  = sigmaSq*(chol(X'*X)\inv(chol(X'*X))');
    beta    = mvnrnd(mu,VarCov,1000);

for mc = 1:1000
     %Get random correlated beta value from normal distribution
    BB = beta(mc,:);
     %Get distribution of swe
    tempSWE.(glacier) = repmat(BB(1),options.mapsize(g,:));

        fields = fieldnames(topo_full.(glacier));
    for f = 1:length(fields)
        param = fields{f};
        val = topo_full.(glacier).(param)*BB(f+1);
        tempSWE.(glacier) = tempSWE.(glacier) + val; 
    end
        tempSWE.(glacier)(tempSWE.(glacier)<0) = 0;

LRtemp.(den)(mc).(glacier) = tempSWE.(glacier);
        
 %Winter balance
fullQbeta.(den).(glacier)(mc) = nanmean(tempSWE.(glacier)(:));

end
    end
end
clock
e = (cputime-t)/60/60

%save('WSMBDistribution.mat','var*','LR*','options','-v7.3')

%% WSMB - Beta & SWE Var

format shortg
clock
t = cputime;
%load Full.mat; load TopoSWE.mat
%clear fullQbetazz

for d = 1:8
    den = options.DenOpt{d};
[ inputSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled);
end

varSIG = [0.027, 0.035, 0.04];
for g = 1:3;    glacier = options.glacier{g};
    varPD.(glacier)   = makedist('Normal','mu',0,'sigma',varSIG(g));
    varSWE.(glacier)  = random(varPD.(glacier), length(inputSWE.S1.(glacier)),1000);   
end

for d = 1%:8
    den = options.DenOpt{d};
    display(den)
for vc = 1%:1000
    for g = 1:3;        glacier = options.glacier{g};
        dataSWE.(den).(glacier)(:,1) = inputSWE.(den).(glacier)(:,1) + varSWE.(glacier)(:,vc);
        I = (dataSWE.(den).(glacier)(:,1)<0);
        dataSWE.(den).(glacier)(I,1) = 0;
    end

% Linear regression

    LRtempbetazz.(den)(vc)  =  LinearRegression( dataSWE.(den), TOPOdata, topo_full );
    
    for g = 1:3;        glacier = options.glacier{g};
        mu      = LRtempbetazz.(den)(vc).coeff{[8,1:7],g};
    
 %Get Beta dist
    X       = struct2table(TOPOdata.(glacier)); X = X{:,:}; X = [ones(length(X),1) X];
    Y       = dataSWE.(den).(glacier)(:,1);
    
    %VarCov matrix from betas
    sigmaSq = sum((Y-X*mu).^2)/(size(X,1)-size(X,2));
    VarCov  = sigmaSq*(chol(X'*X)\inv(chol(X'*X))');
    beta    = mvnrnd(mu,VarCov,1000);

for mc = 1:1000
     %Get random correlated beta value from normal distribution
    BB = beta(mc,:);
     %Get distribution of swe
    tempSWE.(glacier) = repmat(BB(1),options.mapsize(g,:));

        fields = fieldnames(topo_full.(glacier));
    for f = 1:length(fields)
        param = fields{f};
        val = topo_full.(glacier).(param)*BB(f+1);
        tempSWE.(glacier) = tempSWE.(glacier) + val; 
    end
        tempSWE.(glacier)(tempSWE.(glacier)<0) = 0;

 %Winter balance
fullQbetazz.(den).(glacier)(mc,vc) = nanmean(tempSWE.(glacier)(:));

end
    end
end
clock
e = (cputime-t)/60/60
end

%% PLOT -> LR

 %all Gs, one density
    x = 0.2:0.001:0.8;
for o = 1:3
    if      o == 1; data = varB.LR.zz;     t = '\sigma_{SWE} Variability';      f = 'zz';
    elseif  o == 2; data = varB.LR.interp;   t = '\sigma_{\beta} Variability';    f = 'beta';
    elseif  o == 3; data = varB.LR.zzinterp; t = '\sigma_{SWE} and \sigma_{\beta} Variability'; f = 'betaNzz';
    end
figure(1); clf; p = 1;
for d = 1:8;    den = options.DenOpt{d};
for g = 1:3;    glacier = options.glacier{g};
    
        ProbDen = fitdist(data.(den).(glacier)(:),'Normal');
        y       = pdf(ProbDen,x);   
subplot(4,2,p)    
fill(x,y,options.RGB(g,:),'FaceAlpha',0.85, 'EdgeColor', 'none'); hold on
    ylabel('Density'); xlabel('WSMB (m w.e.)');
    xlim([min(x) max(x)])
end
    if p ==8; legend(options.glacier, 'Location','northwest'); end
    title(den);
p = p+1;
end
saveFIG(['WSMB_Distribution',f],12)
end

    
 %With ZZ var vs Without
  den = options.DenOpt{1};
    c = [119, 24, 15; 71, 72, 137; 214, 164, 47]/255;
 figure(4); clf
 bins    = 200;%round(sqrt(length(Qbetazz.(den).G4(:))));
 edges   = linspace(0.3,0.8,bins); 
for g = 1:3
    glacier = options.glacier{g};
subplot(1,3,g)
    histogram(varBzz.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.6, 'FaceColor',c(1,:)); hold on
    histogram(varBbeta.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.6, 'FaceColor',c(2,:)); hold on
    histogram(varBbetazz.(den).(glacier), edges, 'Normalization','probability', 'EdgeColor','none','FaceAlpha',0.6, 'FaceColor',c(3,:)); hold on
   
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); title([glacier,' WSMB Distribution - ',den])
    legend('\sigma_{ZZ}','\beta','\beta and \sigma_{ZZ}','location','best') 
end
    saveFIG('WSMB_compareBetaAndZZBeta')
    
    
 %all Density, one G
    c = cbrewer('qual','Set1',8);
    x = 0.2:0.001:0.8;
for o = 1:3
    if      o == 1; data = varB.LR.zz;     t = '\sigma_{SWE} Variability';      f = 'zz';
    elseif  o == 2; data = varB.LR.interp;   t = '\sigma_{\beta} Variability';    f = 'beta';
    elseif  o == 3; data = varB.LR.zzinterp; t = '\sigma_{SWE} and \sigma_{\beta} Variability'; f = 'betaNzz';
    end
figure(1); clf; 
for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};

    ProbDen = fitdist(data.(den).(glacier)(:),'Normal');
    y       = pdf(ProbDen,x);   
subplot(1,3,g)
fill(x,y,c(d,:),'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on
    ylabel('Density'); xlabel('WSMB (m w.e.)'); title(glacier)
    xlim([min(x) max(x)])
    title(options.glacier{g})
end
    if g ==3; legend(options.DenOpt); end
end
    saveFIG(['WSMB_allDensity_',f],13)    
end

 %all Density, one G
figure(2); clf
c = [147, 148, 150; 37, 37, 38]/255;   
p = 1;
for d = 1:8
    den = options.DenOpt{d};
for g = 1:3
glacier = options.glacier{g};

subplot(8,3,p)
    histogram(fullQbeta.(den).(glacier), 'Normalization','probability',...
        'FaceColor',c(2,:),'EdgeColor','none'); hold on
    histogram(fullQbetazz.(den).(glacier), 'Normalization','probability',...
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
Tbeta.(glacier)(:,:,d)  = varBbeta.(den).(glacier);
Tzz.(glacier)(:,:,d)    = varBzz.(den).(glacier);
Tdlr.(glacier)(:,:,d)     = nanmean(fullLR.(den).(glacier)(:));
Tdsk.(glacier)(:,:,d)     = nanmean(fullSK.(den).(glacier)(:));
Tdbetazz.(glacier)(:,:,d)       = varBbetazz.(den).(glacier);
Tdzzkriging.(glacier)(:,:,d)    = varBkriging.(den).(glacier);
end
end


%dNbetaNzz
figure(1); clf
for g = 3:-1:1
glacier = options.glacier{g};
histogram(Tdbetazz.(glacier)(:),'Normalization','probability',...
        'EdgeColor','none','FaceAlpha',0.75, 'FaceColor',options.RGB(g,:)); hold on
end
ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
legend(options.glacier)
       % saveFIG('WSMB_dNbetaNzz')
       
%dNbetaNzz
figure(2); clf
for g = 3:-1:1
glacier = options.glacier{g};
% histogram(Tdzzkriging.(glacier)(:),'Normalization','probability',...
%         'EdgeColor','none','FaceAlpha',0.75, 'FaceColor',options.RGB(g,:)); hold on
histfit(Tdzzkriging.(glacier)(:)); hold on

end
ylabel('Probability'); xlabel('WSMB (m w.e.)'); 
legend(options.glacier)
       % saveFIG('WSMB_dNbetaNzz')
       

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% KRIGING
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%clear Q
%load TopoSWE.mat
for d = 1:8
    den = options.DenOpt{d};
[ inputSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled);
end

varSIG = [0.027, 0.035, 0.04];
for g = 1:3;    glacier = options.glacier{g};
    varPD.(glacier)   = makedist('Normal','mu',0,'sigma',varSIG(g));
    varSWE.(glacier)  = random(varPD.(glacier), length(inputSWE.S1.(glacier)),1000);   
end

for d = 1:8
    den = options.DenOpt{d};
    display(den)
for vc = 1:100%0
    display(num2str(vc))
    for g = 1:3;        glacier = options.glacier{g};
        dataSWE.(den).(glacier)      = inputSWE.(den).(glacier);
        dataSWE.(den).(glacier)(:,1) = inputSWE.(den).(glacier)(:,1) + varSWE.(glacier)(:,vc);
        I = (dataSWE.(den).(glacier)(:,1)<0);
        dataSWE.(den).(glacier)(I,1) = 0;
    end
    
% Simple kriging

  tempSKzz.(den)(vc) =  KrigingR_G( dataSWE.(den) );

 %Winter balance
for g = 1:3
        glacier = options.glacier{g};      
        Qkrigingzz.(den).(glacier)(vc) = nanmean(tempSKzz.(den)(vc).(glacier).pred(:));
end
end
end

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

save('Kriging1.mat','Qkrigingzz','tempSKzz','-v7.3')
%% PLOT -> Kriging
 %all Gs, one density
bins    = 200;%round(sqrt(length(Qbeta.(den).G4(:))));
edges   = linspace(0.15,0.7,bins); 

data = varBkriging;     t = '\sigma_{ZZ} Variability';      f = 'zz';
figure(1); clf; p = 1;
for d = 1:8
    den = options.DenOpt{d};
for g = 1:3
    glacier = options.glacier{g};
subplot(4,2,p)    
    h(g) = histogram(data.(den).(glacier)(:), edges, 'Normalization','probability',...
        'EdgeColor','none','FaceAlpha',0.5, 'FaceColor',options.RGB(g,:)); hold on
    ylabel('Probability'); xlabel('Winter surface mass balance (m w.e.)');  
end
    legend(h,options.glacier,'Location','northwest')
    title([t,den]);
p = p+1;
end
saveFIG(['WSMB_Distribution_Kriging_separateD_',f],12)

    
 %all Density, one G
    c = cbrewer('qual','Set1',8);
    x = 0:0.001:1.2;
for o = 1:3
    if      o == 1; data = varB.SK.zz;       t = '\sigma_{SWE} Variability';      f = 'zz';
    elseif  o == 2; data = varB.SK.interp;   t = '\sigma_{\beta} Variability';    f = 'beta';
    elseif  o == 3; data = varB.SK.zzinterp; t = '\sigma_{SWE} and \sigma_{\beta} Variability'; f = 'betaNzz';
    end
figure(1); clf; 
for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};

    ProbDen = fitdist(data.(den).(glacier)(:),'Normal');
    y       = pdf(ProbDen,x);   
subplot(1,3,g)
fill([0 x],[0 y],c(d,:),'FaceAlpha',0.5, 'EdgeColor', 'none'); hold on
    ylabel('Density'); xlabel('WSMB (m w.e.)'); title(glacier)
    xlim([min(x) max(x)])
    title(options.glacier{g})
end
    if g ==3; legend(options.DenOpt); end
end
    saveFIG(['WSMB_SK_allDensity_',f],13)    
end    
    
    x = 0:0.001:1.2;
for o = 1:3
    if      o == 1; data = varB.SK.zz;       t = '\sigma_{SWE} Variability';      f = 'zz';
    elseif  o == 2; data = varB.SK.interp;   t = '\sigma_{\beta} Variability';    f = 'beta';
    elseif  o == 3; data = varB.SK.zzinterp; t = '\sigma_{SWE} and \sigma_{\beta} Variability'; f = 'betaNzz';
    end
figure(1); clf; p = 1;
for d = 1:8;    den = options.DenOpt{d};
for g = 1:3;    glacier = options.glacier{g};
    
        ProbDen = fitdist(data.(den).(glacier)(:),'Normal');
        y       = pdf(ProbDen,x);   
subplot(4,2,p)    
fill([0 x],[0 y],options.RGB(g,:),'FaceAlpha',0.85, 'EdgeColor', 'none'); hold on
    ylabel('Density'); xlabel('WSMB (m w.e.)');
    xlim([min(x) max(x)])
end
    if p ==8; legend(options.glacier, 'Location','northwest'); end
    title(den);
p = p+1;
end
saveFIG(['WSMB_Distribution_SK_',f],12)
end

%% Heatmap, spatial variability -> SWE var 
%load WSMBDistribution.mat

%data = varMAP.LR.zzinterp;  method = 'LR';  t = 'betazz';
data = varMAP.SK.zzinterp; method = 'SK'; t = 'krigingzz';

for d = 1:8; den = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
s = size(data.(den)(1).(glacier).pred);

clear F G A H U W X 
n = 1; p = 1;
runs = 100;
for i = 1:runs
    F(n:n+s(1)-1,:) = data.(den)(i).(glacier).pred;
    G(:,p:p+s(2)-1) = data.(den)(i).(glacier).pred;
    n = n+s(1);
    p = p+s(2);
end

F = repmat(F,1,runs);
G = repmat(G,runs,1);

H = F-G;
H = tril(H);
H(H<=0) = NaN;
U = isnan(H);

n=1;
for i = 1:s(1):size(H,1)
for j = 1:s(2):size(H,2)
    A(:,:,n) = H(i:i+s(1)-1,j:j+s(2)-1);
    W(:,:,n) = U(i:i+s(1)-1,j:j+s(2)-1);
    n=n+1;
end
end

for i = size(A,3):-1:1
   X(i) = all(all(W(:,:,i)));
end
A(:,:,X) = [];
    
D.(method).(den).(glacier) = nansum(A,3);
    D.(method).(den).(glacier) = D.(method).(den).(glacier)/max(D.(method).(den).(glacier)(:));
    D.(method).(den).(glacier)(options.mapNaN.(glacier)) = NaN;
end
end

% % SWE Var Map
% den = 'S1';
% for g = 1:3; glacier = options.glacier{g};
%     DD.(glacier) = D_LR.(den).(glacier)/(max(D_LR.(den).(glacier)(:))*0.55);
% end
% 
% figure(1); 
% PlotTopoParameter(DD,'hot','Variability',SWE,'none','nomassB')
%     saveFIG(['SpatialVariabilityMap_SWEVAR_',t,den])

%% Heatmap, spatial variability -> DENSITY var (mean swe)
%load WSMBDistribution.mat

data = varLRbetazz;  t = 'betazz';

clear P K C

for d = 1:8; den = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
    clear P
    for i=1:length(data.(den))
    P(:,:,i) = data.(den).(glacier);
    end
    K.(den).(glacier) = nanmean(P,3);
end
end

for g = 1:3; glacier = options.glacier{g};
s = size(K.(den)(1).(glacier));
clear F G A H U W X 

n = 1; p = 1;
for d = 1:8; den = options.DenOpt{d};
    F(n:n+s(1)-1,:) = K.(den).(glacier);
    G(:,p:p+s(2)-1) = K.(den).(glacier);
    n = n+s(1);
    p = p+s(2);
end

F = repmat(F,1,8);
G = repmat(G,8,1);

H = F-G;
H = tril(H);
H(H<=0) = NaN;
U = isnan(H);

n=1;
for i = 1:s(1):size(H,1)
for j = 1:s(2):size(H,2)
    A(:,:,n) = H(i:i+s(1)-1,j:j+s(2)-1);
    W(:,:,n) = U(i:i+s(1)-1,j:j+s(2)-1);
    n=n+1;
end
end

for i = size(A,3):-1:1
   X(i) = all(all(W(:,:,i)));
end
A(:,:,X) = [];
    
C.(glacier) = nansum(A,3);

end


% DENSITY Var Map
for g = 1:3; glacier = options.glacier{g};
    CC.(glacier) = C.(glacier)/(max(C.(glacier)(:))*0.55);
    CC.(glacier)(options.mapNaN.(glacier)) = NaN;
end

figure(1); 
PlotTopoParameter(CC,'hot','Variability',SWE,'none','nomassB')
    saveFIG(['SpatialVariabilityMap_DENSITY_',t])


%% STD of distributions

for o = 1:3
    if      o == 1; data = varBzz;     Pmax = 53;   t = '\sigma_{ZZ} Variability';       f = 'zz';
    elseif  o == 2; data = varBbeta;   Pmax = 21;   t = '\beta Variability';             f = 'beta';
    elseif  o == 3; data = varBbetazz; Pmax = 53;   t = '\beta and \sigma_{ZZ} Variability'; f = 'betaNzz';
    end

for g = 1:3;    glacier = options.glacier{g};  
for d = 1:8;    den = options.DenOpt{d};
    
    ProbDenLR.(f).(den).(glacier) = fitdist(data.(den).(glacier)(:),'Normal');
sigmaLR.(f)(d,g) = ProbDenLR.(f).(den).(glacier).sigma;
    
end
end
end
    mean(sigmaLR.zz)
    mean(sigmaLR.beta)
    mean(sigmaLR.betaNzz)

for g = 1:3;    glacier = options.glacier{g};  
for d = 1:8;    den = options.DenOpt{d};
    
    ProbDenSK.(f).(den).(glacier) = fitdist(varBsk.(den).(glacier)(:),'Normal');
sigmaSK.(f)(d,g) = ProbDenSK.(f).(den).(glacier).sigma;
    
end
end

%% Ablation area only

%%% MAKE NAN MASKS
 load Basic.mat 
nanMask.G4 = mapNaN.G4;
    nanMask.G4(1:20,:) = 1;
    nanMask.G4(:,1:15) = 1;
    Iy = 40:-2:18;   n = 1;
    for i = 19:30
       nanMask.G4(1:Iy(n),i) = 1; 
       n = n+1;
    end
    Iy = 20:2:34;   n = 1;
    for i = 40:47
       nanMask.G4(1:Iy(n),i) = 1; 
       n = n+1;
    end

nanMask.G2 = mapNaN.G2;
    nanMask.G2(:,100:127) = 1;
    nanMask.G2(95:end,:) = 1;
    Ix = 85:100;   n = 1;
    for i = 100:-1:85
       nanMask.G2(i,Ix(n):end) = 1; 
       n = n+1;
    end
    
nanMask.G13 = mapNaN.G13;
    nanMask.G13(:,95:end) = 1;
    nanMask.G13(110:end,:) = 1;
    nanMask.G13(93:end,1:60) = 1;
    Iy = 50:110;   n = 1;
    for i = 20:80
       nanMask.G13(Iy(n):end,i) = 1; 
       n = n+1;
    end
    
%%% APPLY MASKS TO VAR MAPS
load varWSMB.mat varMAP options
clear ablB
% LR
for u = 1:3
    if      u == 1; uncert = 'zz';
    elseif  u == 2; uncert = 'interp';
    elseif  u == 3; uncert = 'zzinterp';
    end
for d = 1:8; density = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
    
for i = 1:1000
    map = varMAP.(method).(uncert).(density)(i).(glacier);
    map(nanMask.(glacier)) = NaN;
ablB.(method).(uncert).(density).(glacier)(i) = nanmean(map(:));
end

end
end
end

% SK
for u = [1,3]
    if      u == 1; uncert = 'zz';
    elseif  u == 2; uncert = 'interp';
    elseif  u == 3; uncert = 'zzinterp';
    end
for d = 1:8; density = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
    
for i = 1:1000
    map = varMAP.(method).(uncert).(density)(i).(glacier).pred;
    map(nanMask.(glacier)) = NaN;
ablB.(method).(uncert).(density).(glacier)(i) = nanmean(map(:));
end

end
end
end

%%% APPLY MASK TO FULL LR AND SK MAPS

load Full.mat
% LR & SK
for d = 1:8;    den = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};    
    map = fullLR.(den).(glacier);
    map(nanMask.(glacier)) = NaN;
ablB.LR(d,g) = nanmean(map(:));
    map = fullSK.(den).(glacier).pred;
    map(nanMask.(glacier)) = NaN;
ablB.SK(d,g) = nanmean(map(:));
end
end

D = [mean(ablB.LR); mean(ablB.SK)];
bar(D'); title('Ablation only winter balance'); legend('LR','SK')
%%
%%% PLOT

figure(9); clf;
x = 0.15:0.001:0.8;
for o = 4:5
    if      o == 1; data = ablB.LR.zz;      t = '\sigma_{SWE}';     
    elseif  o == 2; data = ablB.LR.interp;  t = '\sigma_{INT}'; 
    elseif  o == 3; data = ablB.LR.zzinterp;  t = '\sigma_{ALL}'; 
        
    elseif  o == 4; data = ablB.SK.zz;      t = '\sigma_{SWE}';     
    elseif  o == 5; data = ablB.SK.zzinterp;  t = '\sigma_{INT}';    
    end
for g = 1:3;     glacier = options.glacier{g};
for d = 1:8;     den = options.DenOpt{d};
    
ProbDen.(den).(glacier) = fitdist(data.(den).(glacier)(:),'Normal');
    y = pdf(ProbDen.(den).(glacier),x);  
    
if  o == 5;    y = [0 y]; x = [x 0];    end

subplot(1,3,o-3) 
fill(x,y,options.RGB(g,:),'FaceAlpha',0.2, 'EdgeColor', 'none'); hold on
%     if      o == 1;         ylabel('LR Density');
%     elseif  o == 4;         ylabel('SK Density');  
%     end
%     if  o == 4||o == 5; xlabel('WSMB (m w.e.)');  end
    
end
end
    xlim([min(x) max(x)])
    title(t);
end