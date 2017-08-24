    %% Load data and calculate
load TopoBMS_MLR.mat SWE options
load Full.mat fullLR fullSK
load Subset.mat

subsets = fieldnames(subsetLR); subsets = sort(subsets); 
    %subsets = subsets([3,1,2,4,5,8,6,7,9,10]);
    T = zeros(8,2); rmseLR = zeros(3,8); sweLR = zeros(3,8);
for g = 1:3;               glacier = options.glacier{g};
for s = 1:length(subsets); sub = subsets{s};
for nn = 1:length(subsetLR)
for d = 1:8;               den = options.DenOpt{d};
    resLR(d) = SampledCell( fullLR.(den) );
    rmseLR(g,d) = sqrt(mean((SWE(g).swe-resLR(d).(glacier)).^2));
    sweLR(g,d)  = nanmean(fullLR.(den).(glacier)(:));

    T(d,1)= nanmean(subsetRmseLR(nn).(sub).(den).(glacier)(:));
    T(d,2)= nanmean(subsetLR(nn).(sub).(den).(glacier)(:));
    
end
    Trmse.(glacier)(nn,s)= mean(T(:,1));
    Tswe.(glacier)(nn,s)= mean(T(:,2));
end
end
end

Nsize = size(Trmse.G4);
    Nhalf = Nsize(1,2)/2;
n = repmat([10:10:100]',1,Nsize(1,2)); 
    n(:,strcmp('centreline',subsets)) = 10:5:55;
    n(:,strcmp('Acentreline',subsets)) = 10:5:55;

rmseLR = mean(rmseLR,2);
sweLR = mean(sweLR,2);    
WBerr = [0.029,0.049,0.032]*1.96;

    %Theta 
for g = 1:3;    glacier = options.glacier{g};
for c = 1:length(subsetSK)
    for f = 1:length(subsets)
    for d = 1:8
        den = options.DenOpt{d};
    theta.(glacier)(c,f,d) = subsetSK(c).(subsets{f}).(den).Model(g).theta;
    end
    end
end
end

for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};
fullTheta(d,g) = fullSK.(den).Model(g).theta;
end
theta.(glacier) = mean(theta.(glacier),3);
end


%% Plot -> RMSE and Winter balance and theta (range length)
    %Accumulation subsets are solid lines, no accumulation is dashed lines
figure(1); clf
    cols = [84 13 110; 238 66 102; 
            255 210 63; 72 226 190; 
            18 96 181; 23 93 8]/255;

for g = 1:3;     glacier = options.glacier{g};
    %RMSE
 subplot(3,3,g)
 for l = 1:Nhalf
    p(l) = plot(n(:,l),Trmse.(glacier)(:,l),'Color',cols(l,:),'LineWidth',3); hold on
    plot(n(:,l+Nhalf),Trmse.(glacier)(:,l+Nhalf),'--','Color',cols(l,:),'LineWidth',3)
 end

     plot([0, max(n(:))],[rmseLR(g),rmseLR(g)],'k--')    
        ylabel('RMSE (m w.e.)'); 
        title(glacier)
        if      g==1; ylim([0.14 0.24])
        elseif  g==2; ylim([0.10 0.17]);
        elseif  g==3; ylim([0.08 0.14]);     
            legend(p,subsets(1:6));
        end
        set(gca,'YTick',(0:0.01:1))
        
    %Winter balance    
 subplot(3,3,g+3)
 for l = 1:Nhalf
    plot(n(:,l),Tswe.(glacier)(:,l),'Color',cols(l,:),'LineWidth',3); hold on
    plot(n(:,l+Nhalf),Tswe.(glacier)(:,l+Nhalf),'--','Color',cols(l,:),'LineWidth',3)
 end
    %plot([0, max(n(:))],[sweLR(g), sweLR(g)],'k--')   
    h = fill([0, max(n(:)), max(n(:)),0],...
         [sweLR(g)+WBerr(g), sweLR(g)+WBerr(g), sweLR(g)-WBerr(g), sweLR(g)-WBerr(g)],...
         [161, 162, 163]/255);  set(h,'facealpha',.5); set(h,'EdgeColor','none')
     

        ylabel('LR Winter balance (m w.e.)'); %xlabel('Sample Size')
        ylim([0.2 0.8]);
        
    %Range length (theta)    

subplot(3,3,g+6)
for l = 1:Nhalf
    plot(n(:,l),theta.(glacier)(:,l),'Color',cols(l,:),'LineWidth',3); hold on
    plot(n(:,l),theta.(glacier)(:,l+Nhalf),'--','Color',cols(l,:),'LineWidth',3); hold on
end
    ylabel('\theta (m)'); xlabel('Sample size')
%fill([min(T) min(T)],[max(T) max(T)],[207, 207, 209]);
    
    %plot([0 100],[T(g) T(g)],'k', 'LineWidth',3)
        
%     insetx = [.27, .55, .83];
%     axes('Position',[insetx(g) .67 .06 .1]); box on
%     histogram(theta.(glacier)(:,:),15, 'FaceColor', options.RGB(g,:))
%         ylabel('Frequency'); xlabel('\theta (m)')

end

    saveFIG('SubsetInterpSizeCompile',14)

% Display subset with highest RMSE    
for g = 1:3;     glacier = options.glacier{g};
        L = max(mean(Trmse.(glacier)(4:10,:))) == mean(Trmse.(glacier)(4:10,:));
    display([glacier, 'max: ',subsets{L}])
end
% Display subset with lowest RMSE
for g = 1:3;     glacier = options.glacier{g};
        L = min(mean(Trmse.(glacier)(4:10,:))) == mean(Trmse.(glacier)(4:10,:));
    display([glacier, 'min: ',subsets{L}])
end



%% Variogram
load Variogram.mat
load TopoSWE.mat options

    N = zeros(3,1); S = zeros(3,1); R = zeros(3,1);
for g = 1:3;     glacier = options.glacier{g};
    N(g) = mean(nugget.(glacier)(:));
    R(g) = nanmean(range.(glacier)(:));
    S(g) = mean(sill.(glacier)(:));
    
% subplot(3,3,g)
%     histogram(nugget.(glacier)(:)); ylabel('Nugget')
% subplot(3,3,g+3)
%     histogram(range.(glacier)(:)); ylabel('Range')
% subplot(3,3,g+6)
%     histogram(sill.(glacier)(:)); ylabel('Sill')
   
    x = 0:round(R(g)*1.3,-3);
    y = N(g) + ( S(g)*( 1.5*(x/R(g)) - 0.5*(x/R(g)).^3).*(x <= R(g)) + S(g)*(x>R(g)));
figure(2);
subplot(1,3,g)
    plot(x,y)
    title(glacier); ylabel('Semivariance'); xlabel('Distance (m)')
end

%%

