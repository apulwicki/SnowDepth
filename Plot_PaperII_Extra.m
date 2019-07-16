figure(1); clf

% title_beta = {'Elevation', 'Aspect', 'Slope', 'Northness', 'Curvature', 'Sx'};
title_beta = {'Elevation', 'Slope', 'Curvature', 'Sx'};

num_models = 200;

%McGrath Values
% elev_range  = [0.585, 0.810];
% dc_range    = [0, 0.015];
% aspect_range = [0, 0.093];
% slope_range = [-0.09, 0.277];
% N_range     = [-0.414, 0.172];
% curv_range  = [-0.077, 0.02];
% sx_range    = [-0.294, 0.260];
% Our values
elev_range  = [0.005, 0.118];
dc_range    = [0, 0.015];
aspect_range = [0, 0.093];
slope_range = [0, 0.038];
N_range     = [-0.414, 0.172];
curv_range  = [-0.028, 0.035];
sx_range    = [-0.055, 0.04];

% norm_dist = [mean(elev_range), std(elev_range);
%             mean(aspect_range), std(aspect_range);
%             mean(slope_range), std(slope_range);
%             mean(N_range), std(N_range);
%             mean(curv_range), std(curv_range);
%             mean(sx_range), std(sx_range)];
norm_dist = [mean(elev_range), std(elev_range);
            mean(slope_range), std(slope_range);
            mean(curv_range), std(curv_range);
            mean(sx_range), std(sx_range)];

for i=1:4
   subplot(2,2,i)
   coeffs = normrnd(norm_dist(i,1),norm_dist(i,2), [1,num_models]);
   histogram(coeffs,100,'EdgeColor','none','Normalization','pdf'); hold on
   
   x = [norm_dist(i,1)-3*norm_dist(i,2):0.01:norm_dist(i,1)+3*norm_dist(i,2)];
   norm = normpdf(x,norm_dist(i,1),norm_dist(i,2));
   plot(x,norm)
   
   title(title_beta(i))
   xlabel('\beta')
   ylabel('Frequency')
   xlim([min(x) max(x)])
end
	saveFIG_HP('Syn_BetaDist_4var_PulwickiCoeff',2,17)


%%

    load PaperII_AblationArea.mat AblationArea
    load PaperII_RegularSampling.mat 
    load Full.mat fullLR
    run OPTIONS

% Blank out accum area
for g = 1:3;    glacier = options.glacier{g};
    abla_area.(glacier) = fullLR.S2.(glacier);
    abla_area.(glacier)(AblationArea.(glacier)==-0.1)=-0.1;
end

figure(1); clf
[ha, pos] = tight_subplot(5,6,[0.04 0],[0 0.05],[0 0]);
s=1;
for n=1:10
for g=1:3; glacier=options.glacier{g};
axes(ha(s));    
h(s) = imagesc(snowdist_model(n).(glacier)); hold on
    set(h(s),'alphadata',~isnan(snowdist_model(n).(glacier)));
    axis square; axis off;
    caxis([0 1.5])
    if g==2; title(['Model ',num2str(n)]); end
    s=s+1;
end
end
	saveFIG_HP('Syn_SnowDist',2,17)

    
    %%
    load PII_SynDistributions

for m=0:19
figure(1); clf
s=1;
for n=1:10
for g=1:3; glacier=options.glacier{g};
    subplot(5,6,s)
    model_num = n+10*m;
    snow_dist = snowdist_model(model_num).(glacier);
    snow_dist(AblationArea.(glacier)==-0.1)=-0.1;
    h = imagesc(snow_dist); hold on
    
    Cmap = parula;      Cmap(1,:) = [199, 201, 204]/255;
    colormap(Cmap);     set(h,'alphadata',~isnan(snow_dist));
    axis square; axis off;
    caxis([-0.1 1])
    if g==2; title(['Model ',num2str(model_num)]); end 
    s=s+1;
end
end
pause
end

%% Calculate Bw in ablation area only for scaling
load Full.mat fullLR
load PaperII_AblationArea.mat AblationArea

for g=1:3; glacier=options.glacier{g};

    fullLR.S2.(glacier)(AblationArea.(glacier)==-0.1)=NaN;
%     imagesc(fullLR.S2.(glacier))
    abl_meanBw(g) = nanmean(fullLR.S2.(glacier)(:));
    
end

%% RMSE of diff scenarios

% clear
load PaperII_Syntheic.mat
load PaperII_RegularSampling
load PaperII_AblationArea.mat AblationArea AccumulationArea
% load PaperII_FinalLRruns.mat snowdist_model coeffSYN

file_path = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/';
% file_path = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/';

sampling_number = 200;
num_models = 1000;
nn = 6:45;

% Get synthetic snow distributions

Bw_abl = [0.59, 0.34, 0.27]; % ablation area only
Bw_full = [0.59, 0.57, 0.38]; % whole glacier
Bw_obs = Bw_full;
    p_acc = [0.36,    0.3,        0.308];
    p_abl = [0.3422,  0.353,     0.3692];

% Beta coeff ranges
% McGrath values
elev_range  = [0.585, 0.810];
dc_range    = [0, 0.015];
aspect_range = [0, 0.093];
slope_range = [-0.09, 0.277];
N_range     = [-0.414, 0.172];
curv_range  = [-0.077, 0.02];
sx_range    = [-0.294, 0.260];

% Our values
% elev_range  = [0.005, 0.118];
% dc_range    = [0, 0.015];
% aspect_range = [0, 0.093];
% slope_range = [0, 0.038];
% N_range     = [-0.414, 0.172];
% curv_range  = [-0.028, 0.035];
% sx_range    = [-0.055, 0.04];


intercept_range   = [0.6205, 0.2627, 0.2354];

% Normal random beta values
intercept   = normrnd(mean(intercept_range), std(intercept_range), [num_models,3]);
elev_beta   = normrnd(mean(elev_range), std(elev_range), [1,num_models]);
dc_beta     = normrnd(mean(dc_range), std(dc_range), [1,num_models]);
aspect_beta = normrnd(mean(aspect_range), std(aspect_range), [1,num_models]);
slope_beta  = normrnd(mean(slope_range), std(slope_range), [1,num_models]);
N_beta      = normrnd(mean(N_range),std(N_range),[1,num_models]);
curv_beta   = normrnd(mean(curv_range), std(curv_range), [1,num_models]);
sx_beta     = normrnd(mean(sx_range),   std(sx_range),   [1,num_models]);

coeffSYN = [elev_beta', slope_beta', curv_beta', sx_beta'];


    % Remove dc, aspect and Northness
    for g = 1:3;    glacier = options.glacier{g};
    topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
    end

for m = 1:num_models

for g = 1:3;    glacier = options.glacier{g};
    syn_model   = repmat(intercept(m,g), options.mapsize(g,:));
    betaCoeff   = coeffSYN(m,:);   
    topoCoeff   = fieldnames(topo_full.G4);
    
    for c = 1:length(betaCoeff)
        param     = topoCoeff{c};
        sweT      = topo_full.(glacier).(param)*betaCoeff(c);
        syn_model = syn_model + sweT;
        
    end
        
    syn_model(syn_model<0) = 0; %Set min to 0
   
    syn_model(AblationArea.(glacier)==-0.1)=NaN;
    scale_to_obs = nanmean(syn_model(:))/Bw_obs(g);

    snowdist_model(m).(glacier) = syn_model/scale_to_obs;    
    
end
end

    % Topo
    for g = 1:3;    glacier = options.glacier{g};
    for c = 1:length(betaCoeff);   param = topoCoeff{c};
    topo_temp.(param).(glacier) = topo_full.(glacier).(param);
    end
    topo.(glacier) = zeros(2,2);
    end

    for g = 1:3;    glacier = options.glacier{g};
        for c = 1:length(betaCoeff);   param = topoCoeff{c};
        tt = SampledCell(topo_temp.(param));
        T(:,c) = tt.(glacier);
        end
        topo.(glacier) = T;
        clear T
    end
    
% All data bw
for m = 1:num_models
    true_swe = SampledCell(snowdist_model(m));
    for g = 1:3;    glacier = options.glacier{g};
        swe	    = true_swe.(glacier);
        Xt      = topo.(glacier);
        X       = [ones(size(Xt,1),1), Xt];
%         X       = X(all(~isnan(X),2));

        coeffs = regress(swe, X);
        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;

            
            % RMSE
        sampledtemp     = swe_pred(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        syn_measure     = diag(sampledtemp);
        
        real_measure_tmp = true_swe.(glacier)(~isnan(syn_measure));
        syn_measure_tmp = syn_measure(~isnan(syn_measure));
        real_measure_tmp2 = true_swe.(glacier)(~isnan(real_measure_tmp));
        syn_measure_tmp2 = syn_measure(~isnan(real_measure_tmp));
               
        rmse_alldata.(glacier)(m) = sqrt(mean((syn_measure_tmp2-real_measure_tmp2).^2));
    end
end


for g = 1:3;    glacier = options.glacier{g};
    display(mean(rmse_alldata.(glacier)))
end

 
%%

mc = 1;
ss=45;

figure(1); clf
n=1;
for g = 1:3;    glacier = options.glacier{g};
 

subplot(3,7,n)
h=imagesc(snowdist_model(mc).(glacier));
set(h,'alphadata',~isnan(snowdist_model(mc).(glacier)));
caxis([0 1]); caxis(caxis);
title(['True, B_w=', num2str(round(nanmean(snowdist_model(mc).(glacier)(:)),2))])
axis off

    for p = 1:length(namesP)    
subplot(3,7,n+p)
data = predLR(ss).(namesP{p})(mc).(glacier);
h=imagesc(data);
set(h,'alphadata',~isnan(data));
title({[namesP{p},', B_w=', num2str(round(mean(synBw_rand.(namesP{p}).(glacier)(ss,mc,:)),2))],...
            ['RMSE=',num2str(round(mean(synRMSE_rand.(namesP{p}).(glacier)(ss,mc,:)),2))]})
caxis([0 1]); caxis(caxis);
axis off 

% display(glacier)
% display(coeffsLR(ss).(namesP{p})(mc).(glacier))
% display(R2.(namesP{p}).(glacier)(ss,mc,sampleN))
% display(synRMSE_rand.(namesP{p}).(glacier)(ss,mc,sampleN))
    end

n=n+7;
end


%%
figure(1); clf
ela_ind = [1 7; 8 15; 16 23];
    ela_m = [0 0; -0.04*10^5 -0.0088*10^6; 0.04*10^5 0.0086*10^6];
for g = 1:3; glacier = options.glacier{g};
subplot(1,3,g)
if g == 3; fill(options.rig.(glacier)(1:304,1),options.rig.(glacier)(1:304,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
else; fill(options.rig.(glacier)(:,1),options.rig.(glacier)(:,2),...
                [199, 201, 204]/255, 'EdgeColor','none'); hold on
end
% plot(options.rig.ELA(ela_ind(g,1):ela_ind(g,2),1)+ela_m(g,1), options.rig.ELA(ela_ind(g,1):ela_ind(g,2),2)+ela_m(g,2),'-k');
plot(option.ELA.(glacier)(:,1),option.ELA.(glacier)(:,2),'-k') 
axis off; axis equal
end

%%

%% Calculate best RMSE using all data


% load PaperII_FinalLRruns
for mc = num_models

%Additing noise to distributed WB
real_measure = SampledCell(snowdist_model(mc));

topoparam = fieldnames(topo_full.G4);
for num_param = 1:length(topoparam)
    for g = 1:3;        glacier = char(options.glacier(g));
    tmp_topo.(glacier) = topo_full.(glacier).(topoparam{num_param});
    end    
    tmp_sampled.(topoparam{num_param}) = SampledCell(tmp_topo);
end
for g = 1:3;        glacier = char(options.glacier(g));
    for num_param = 1:length(topoparam)
sampled_topo.(glacier).(topoparam{num_param}) = tmp_sampled.(topoparam{num_param}).(glacier);
    end
end        
        % BASIC LR        
for g = 1:3;        glacier = char(options.glacier(g));
        display(mc)

    for n=6:length(real_measure.(glacier))
        swe	    = real_measure.(glacier)(:,1);
            swe     = swe + normrnd( 0, options.zzstd(g), size(swe,1),1);
            swe(swe<0) = 0;
        Xt      = struct2array(sampled_topo.(glacier));
        X       = [ones(length(Xt),1), Xt];
        
        X       = X(1:n,:); swe = swe(1:n);
        
        coeffs = regress(swe, X);
        tmp_best_coeffsLR.(glacier)(mc,:)   = coeffs;

        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for num_param = 1:length(betaCoeff)
            param      = topoCoeff{num_param};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_param);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;
        if only_ablation
            swe_pred = swe_pred.*AblationArea.(glacier);
        end        
        
        rmse.(glacier)(n,mc) = sqrt(nanmean(nanmean((snowdist_model(mc).(glacier)-swe_pred).^2)));
    end
    
end   
end

%%
figure(1); clf

for g = 1:3;        glacier = char(options.glacier(g));
   subplot(1,3,g)
%    plot(rmse.(glacier),'Color',[199, 201, 204]/255); hold on
   plot(6:length(real_measure.(glacier)), mean(rmse.(glacier)(6:end,:),2),'k'); hold on
   plot([0,length(real_measure.(glacier))],[best_rmseLR.(glacier) best_rmseLR.(glacier)],'--k')
   ylim([0 0.5])
   xlabel('N')
   ylabel('RMSE (m w.e.)')
   title(glacier)
end

%%
mc = 1;
ss=45;

figure(1); clf
n=0;
namesP = {'Centreline', 'CentreTransect', 'Circle', 'Hourglass', 'HourCircle', 'RandomSafe'};
for g = 1:3;    glacier = options.glacier{g};
 

% subplot(3,7,n)
% h=imagesc(snowdist_model(mc).(glacier));
% set(h,'alphadata',~isnan(snowdist_model(mc).(glacier)));
% caxis([0 1]); caxis(caxis);
% title(['True, B_w=', num2str(round(nanmean(snowdist_model(mc).(glacier)(:)),2))])
% axis off

    for p = 1:length(namesP)    
subplot(3,6,n+p)
data = mean(pred_bw.(namesP{p}).(glacier),3);
h=imagesc(data);
set(h,'alphadata',~isnan(data));
title({[namesP{p},', B_w=', num2str(round(nanmean(data(:)),2))]})%,...
%             ['RMSE=',num2str(round(mean(synRMSE_rand.(namesP{p}).(glacier)(ss,mc,:)),2))]})
caxis([0 1]); caxis(caxis);
axis off 

% display(glacier)
% display(coeffsLR(ss).(namesP{p})(mc).(glacier))
% display(R2.(namesP{p}).(glacier)(ss,mc,sampleN))
% display(synRMSE_rand.(namesP{p}).(glacier)(ss,mc,sampleN))
    end

n=n+6;
end


%%

beta_lab = {'elevation','slope','curvature','Sx'};
nBins = 12;

for n = 1:200
 figure(1); clf
for i =1:4
subplot(2,2,i)
histogram(topo_abl.G13.(beta_lab{i}),nBins,'Normalization','probability'); hold on
histogram(topo_sampled.G13.(beta_lab{i}),nBins,'Normalization','probability'); hold on
histogram(topo_hist.Hourglass(:,i,n),nBins,'Normalization','probability')
histogram(topo_hist.Centreline(:,i,n),nBins,'Normalization','probability')

ylim([0 0.5]); xlim([-10,10])
legend('Full Ablation Area', 'Measured locations','Hourglass','Centreline')
title(beta_lab(i))

end
pause
end