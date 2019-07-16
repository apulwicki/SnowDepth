%% Get synthetic snow distributions
clear
load PaperII_Syntheic.mat topo_full
load PaperII_AblationArea.mat AblationArea AccumulationArea
run OPTIONS

num_models = 200;
only_ablation = true;
save_name = 'PaperII_SynSnowDistModel_Pul_Abl_DDA.mat';


Bw_abl = [0.59, 0.34, 0.27]; % ablation area only
Bw_full = [0.59, 0.57, 0.38]; % whole glacier
if only_ablation
    Bw_obs = Bw_abl;
else
    Bw_obs = Bw_full;
end
    p_acc = [0.36,    0.3,        0.308];
    p_abl = [0.3422,  0.353,     0.3692];


% Beta coeff ranges
% McGrath values
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

% coeffSYN = [elev_beta', dc_beta', aspect_beta', slope_beta', N_beta', curv_beta', sx_beta'];
% coeffSYN = [elev_beta', aspect_beta', slope_beta', N_beta', curv_beta', sx_beta'];
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
    
    syn_model = syn_model + abs(min(syn_model(:)));% + ...
%                 normrnd(0, options.zzstd(g), size(syn_model,1), size(syn_model,2));
    syn_model(syn_model<0) = 0;
    
    if only_ablation
    syn_model(AblationArea.(glacier)==-0.1)=NaN;
    end
    
    syn_model(AccumulationArea.(glacier)==1) = syn_model(AccumulationArea.(glacier)==1)*p_acc(g)./p_abl(g);
    
    scale_to_obs = nanmean(syn_model(:))/Bw_obs(g);

    snowdist_model(m).(glacier) = syn_model/scale_to_obs;

end
end

save(save_name)
% save(save_name,'snowdist_model','coeffSYN')


%% Check that Bw is accurate

clear
load PaperII_SynSnowDistModel_Pul_abl_noise
options.glacier = {'G4','G2','G13'};

Bw = zeros(200,3);
for i=1:200
   for g = 1:3;    glacier = options.glacier{g};
   Bw(i,g) = nanmean(snowdist_model(i).(glacier)(:));
   end
end

%% Plot snow distributions
load PaperII_SynSnowDistModel_Pul_abl_noise

i=40;

figure(1); clf
subplot(1,3,1)
h = imagesc(snowdist_model(i).G4);
set(h,'alphadata',~isnan(snowdist_model(i).G4));
caxis([0 1]); caxis(caxis);
axis off
subplot(1,3,2)
h = imagesc(snowdist_model(i).G2);
set(h,'alphadata',~isnan(snowdist_model(i).G2));
caxis([0 1]); caxis(caxis);
axis off
subplot(1,3,3)
h = imagesc(snowdist_model(i).G13);
set(h,'alphadata',~isnan(snowdist_model(i).G13));
caxis([0 1]); caxis(caxis);
axis off
 
%%
%%
clear 
close all
clc

%load PII_SynDistributions
load PaperII_SynSnowDistModel_McG_accum.mat
%load PaperII_SynSnowDistModel_Pul_accum.mat

options.glacier = {'G4','G2','G13'};

for m=0:19
figure
s=1;
for n=1:10
 for g=1:3; glacier=options.glacier{g};
   subplot(5,6,s)
   model_num = n+10*m;
   snow_dist = snowdist_model(model_num).(glacier);

   h = imagesc(snow_dist); hold on
   set(h,'alphadata',~isnan(snowdist_model(model_num).(glacier)));
   axis square; axis off;
   caxis([0 1.0])
   if g==2; title(['Model ',num2str(model_num)]); end
   s=s+1;

   nanmean(nanmean(snowdist_model(model_num).(glacier)))
   pause

 end
end
end


%%
a = magic(5);
a(1,1) = NaN;

display(nanmean(nanmean(a)))
display(nanmean(a(:)))
