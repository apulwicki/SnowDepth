figure(1); clf

% title_beta = {'Elevation', 'Aspect', 'Slope', 'Northness', 'Curvature', 'Sx'};
title_beta = {'Elevation', 'Slope', 'Curvature', 'Sx'};

num_models = 200;
elev_range  = [0.585, 0.810];
dc_range    = [0, 0.015];
aspect_range = [0, 0.093];
slope_range = [-0.09, 0.277];
N_range     = [-0.414, 0.172];
curv_range  = [-0.077, 0.02];
sx_range    = [-0.294, 0.260];
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
	saveFIG_HP('Syn_BetaDist_4var',2,17)


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