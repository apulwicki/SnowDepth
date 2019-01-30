load PaperII_Syntheic.mat
load PaperII_FinalLRruns.mat snowdist_model coeffSYN intercept

file_path = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/';
%file_path = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/';

num_models = 10;

% Remove dc, aspect and Northness
    for g = 1:3;    glacier = options.glacier{g};
    topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
    end

for mc = 1:num_models

real_measure = SampledCell(snowdist_model(mc));

    %Add some noise
%     WBinputN = WBnoise(WBinput.(namesP{p}),'low');
    
%         % BMA
%         swe_input.(glacier) = WBinputN.(glacier);
%         topo_input = TOPOinput.(namesP{p});
% 
%         cd BMS
%         [BMSinit, BMSres] = BMS_R(swe_input, topo_input);
%         cd ..
% 
%         for g = 3;        glacier = char(options.glacier(g));
%         BMS.(glacier) = BMSinit.(glacier)(:,1);   
%         coeffsBMA(ss).(namesP{p})(mc).(glacier)   = BMS.(glacier){1:8,1};
% 
%         swe_pred = repmat(BMS.(glacier){8,1}, options.mapsize(g,:));
%         betaCoeff = BMS.(glacier){1:7,1};    topoCoeff = fieldnames(topo_full.G4);
%         for num_models = 1:length(betaCoeff)
%             param      = topoCoeff{num_models};
%             sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
%             swe_pred   = swe_pred + sweT;
%         end
%         swe_pred(swe_pred<0) = 0;
%         swe_pred = swe_pred.*AblationArea.(glacier); 
%         predBMA(ss).(namesP{p})(mc).(glacier) = swe_pred;
%         
%         sampledtemp             = swe_pred(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
%         syn_measure.(glacier)   = diag(sampledtemp);
%         
%         real_measure_tmp.(glacier) = real_measure.(glacier)(~isnan(syn_measure.(glacier)));
% 
%             syn_measure_tmp.(glacier) = syn_measure.(glacier)(~isnan(syn_measure.(glacier)));
%         rmseBMA.(namesP{p}).(glacier)(ss,mc) = sqrt(mean((syn_measure_tmp.(glacier)-real_measure_tmp.(glacier)).^2));
%         end
        
        
        % BASIC LR        
        for g = 1:3;        glacier = char(options.glacier(g));

            swe_tmp   = snowdist_model(mc).(glacier)(:);
        swe       = swe_tmp(~isnan(swe_tmp));
        topoCoeff = fieldnames(topo_full.G4); Xt = zeros(length(swe),length(topoCoeff));
        for t=1:length(topoCoeff)
            Xt_temp = topo_full.(glacier).(topoCoeff{t})(:);
            Xt_temp = Xt_temp(~isnan(swe_tmp));
            Xt(:,t) = Xt_temp;
        end
        
        X       = [ones(length(Xt),1), Xt];

        coeffs = regress(swe, X);
        coeffsLR(mc).(glacier)   = coeffs;

        swe_pred = repmat(coeffs(1), options.mapsize(g,:));
        betaCoeff = coeffs(2:end);    topoCoeff = fieldnames(topo_full.G4);
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
        swe_pred(swe_pred<0) = 0;
        swe_pred = swe_pred.*AblationArea.(glacier); 
        predLR(mc).(glacier) = swe_pred;
        
%         sampledtemp     = swe_pred(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
%         syn_measure     = diag(sampledtemp);
%         
%         real_measure_tmp = real_measure.(glacier)(~isnan(syn_measure));
%         syn_measure_tmp = syn_measure(~isnan(syn_measure));
%         
%         rmseLR.(namesP{p}).(glacier)(ss,mc) = sqrt(mean((syn_measure_tmp-real_measure_tmp).^2));
        
        % RMSPE 100%*sum(abs(obs - sim)/obs)
%         x = abs(syn_measure_tmp-real_measure_tmp)./real_measure_tmp;
%         x = x(~isinf(x));
%         x = x(~isnan(x));
%         mpeLR.(namesP{p}).(glacier)(ss,mc) = 100*sum(x)/length(syn_measure_tmp);

        end
        
end


%% Plot coeffs

m = 1;
for i=1:10
    for g = 3;        glacier = char(options.glacier(g));
     subplot(5,2,m)
%     plot_coeffs = [coeffSYN(i,:)',coeffsLR(nn).(pattern)(i).(glacier)(2:end), coeffsBMA(nn).(pattern)(i).(glacier)(1:7);
%                     intercept(g), coeffsLR(nn).(pattern)(i).(glacier)(1), coeffsBMA(nn).(pattern)(i).(glacier)(8)];
    plot_coeffs = [coeffSYN(i,:)',coeffsLR(i).(glacier)(2:end);
                    intercept(g), coeffsLR(i).(glacier)(1)];
   bar(plot_coeffs)
    params={'z','\alpha','m','N','\kappa','Sx','intercept'};
    set(gca,'xticklabel',params)
    ylabel('\beta')
    title(['Model ',num2str(i),' ',glacier])
    m=m+1;
%     legend('Model','Simple LR','BMA')
    end
end

%% 

m = 1;
for i=1:5
    for g = 3;        glacier = char(options.glacier(g));
    subplot(5,2,m)
    actual = snowdist_model(i).(glacier).*AblationArea.(glacier);
    imagesc(actual); colorbar
        caxis([0 2]); caxis(caxis);
    title({['Model ',num2str(i)],['Model ', glacier, ' B_W=', num2str(nanmean(actual(:)))]})
    axis('off')
    
    subplot(5,2,m+1)
    imagesc(predLR(i).(glacier)); colorbar
        caxis([0 2]); caxis(caxis);
    title({['Model ',num2str(i)],['Simple LR ', glacier, ' B_W=', num2str(nanmean(predLR(i).(glacier)(:)))]})
    axis('off')
    
%     subplot(5,3,m+2)
%     imagesc(predBMA(nn).(pattern)(i).(glacier)); colorbar
%         caxis([0 2]); caxis(caxis);
%     title({['Model ',num2str(i)],['BMA ', glacier, ' B_W=', num2str(nanmean(predBMA(nn).(pattern)(i).(glacier)(:)))]})
%     axis('off')
    
    m=m+2;
%     legend('Simple LR','BMA')
    end
end