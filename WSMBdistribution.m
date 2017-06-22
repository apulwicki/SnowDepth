
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
bins    = 200;%round(sqrt(length(Qbeta.(den).G4(:))));
edges   = linspace(0.3,0.8,bins); 
for o = 1%1:3
    if      o == 1; data = varBzz;     t = '\sigma_{ZZ} Variability';      f = 'zz';
    elseif  o == 2; data = varBbeta;   t = '\beta Variability';             f = 'beta';
    elseif  o == 3; data = varBbetazz; t = '\beta and \sigma_{ZZ} Variability'; f = 'betaNzz';
    end
figure(o); clf; p = 1;
for d = 1%:8
    den = options.DenOpt{d};
for g = 1:3
    glacier = options.glacier{g};
subplot(4,2,p)    
    h(g) = histogram(data.(den).(glacier)(:), edges, 'Normalization','probability',...
        'EdgeColor','none','FaceAlpha',0.5, 'FaceColor',options.RGB(g,:)); hold on
    ylabel('Probability'); xlabel('Winter surface mass balance (m w.e.)');  
end
    legend(h,options.glacier)
    title([t,den]);
p = p+1;
end
%saveFIG(['WSMB_Distribution',f],16)
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
for o = 1%:3
    if      o == 1; data = varBzz;     t = '\sigma_{ZZ} Variability';      f = 'zz';
    elseif  o == 2; data = varBbeta;   t = '\beta Variability';             f = 'beta';
    elseif  o == 3; data = varBbetazz; t = '\beta and \sigma_{ZZ} Variability'; f = 'betaNzz';
    end
figure(o); clf; 
for g = 1:3
glacier = options.glacier{g};

for d = 1%:8
    den = options.DenOpt{d};

subplot(1,3,g)
    histogram(data.(den).(glacier), 200,'Normalization','probability',...
        'FaceColor',c(d,:),'EdgeColor','none','FaceAlpha',0.5); hold on
    ylabel('Probability'); xlabel('WSMB (m w.e.)'); title(glacier)
    xlim([0.3 0.8])
    title([t,den])
end
    legend(options.DenOpt)
end
    %saveFIG(['WSMB_allDensityBeta_',f])    
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
       


%% KRIGING

format shortg
clock
t = cputime;

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
clock
e = (cputime-t)/60/60

for g = 1:3; glacier = options.glacier{g};
for i = 1:100
Qupper.(glacier)(i) = nanmean(tempSK.S1(i).(glacier).upper95(:));
Qpred.(glacier)(i) = nanmean(tempSK.S1(i).(glacier).pred(:));
Qlower.(glacier)(i) = nanmean(tempSK.S1(i).(glacier).lower95(:));
end
end

%save('Kriging1.mat','Qkriging','tempSK','-v7.3')
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

%% Heatmap, spatial variability -> SWE var (one density)
%load WSMBDistribution.mat

%data = varLRbetazz;  t = 'betazz';
data = varSKzz;  t = 'krigingzz';

for d = 1; den = options.DenOpt{d};
for g = 1:3; glacier = options.glacier{g};
s = size(data.(den)(1).(glacier));

clear F G A H U W X 
n = 1; p = 1;
runs = 100;
for i = 1:runs
    F(n:n+s(1)-1,:) = data.(den)(i).(glacier);
    G(:,p:p+s(2)-1) = data.(den)(i).(glacier);
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
    
D_SK.(den).(glacier) = nansum(A,3);

end
end

% % SWE Var Map
% den = 'S1';
% for g = 1:3; glacier = options.glacier{g};
%     DD.(glacier) = D_LR.(den).(glacier)/(max(D_LR.(den).(glacier)(:))*0.55);
%     DD.(glacier)(options.mapNaN.(glacier)) = NaN;
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
