type = 'Acentreline';
for g = 1:3;
glacier = options.glacier{g};
Nfull = 10:10:100;
for c = 9:-1:1
    figure(1); clf
    n = Nfull(1:c);
    dd = 7;
den = options.DenOpt{dd};
% Glacier
subplot(3,3,[4,5,7,8])
    h = imagesc(subsetLR(c).(type).(den).(glacier)); hold on
    set(h,'alphadata',~isnan(subsetLR(c).(type).(den).(glacier)))
    col = colorbar; ylabel(col, 'SWE (m w.e.)')
    caxis([0 1.2]); caxis(caxis);
    axis off; axis equal
    
        minE = min(options.rig.(glacier)(:,1));
        minN = min(options.rig.(glacier)(:,2)); 
        Ng = (options.rig.(glacier)(:,2) - minN)/40;  Ng = (max(Ng)-Ng);
                %Ng = Ng + (options.mapsize(3,1)-max(Ng));
        E = (subsetSWE(c).(type).(den).(glacier)(:,2)-minE)/40;
        Na = (subsetSWE(c).(type).(den).(glacier)(:,3)-minN)/40; N = max(Ng)-Na;
    plot(E,N,'k.', 'MarkerSize',5); hold on
   
% Ref Glacier

subplot(3,3,2)
    h = imagesc(fullLR.(den).(glacier)); hold on
    set(h,'alphadata',~isnan(fullLR.(den).(glacier)))
    axis off; axis equal
    caxis([0 1.2]); caxis(caxis);

title(['LR ',type, ' ',glacier])
    
% SWE vs n
    meanLR.(glacier) = nanmean(stackSWE.MLR.(glacier)(:));
    for cc = 1:length(n);
    for d = 1:8 
       den = DenOpt{d};
       stack_swen.(glacier)(d,cc) = nanmean(subsetLR(cc).(type).(den).(glacier)(:));
    end
    end
subplot(3,3,3)
    if c==1
    plot(repmat(n,8,1),stack_swen.(glacier)(:,1:c),'.'); hold on
    else
    plot(n,stack_swen.(glacier)(:,1:c)','LineWidth',2); hold on
    end
    plot([0, max(Nfull)],[meanLR.(glacier), meanLR.(glacier)],'k--')
        xlabel('Sample size'); ylabel('Mean SWE')
        %legend([DenOpt, {'All'}],'Location','best')
        xlim([0, max(Nfull)])
        if c == 9; 
        ymaxswe = ylim;
        else
        ylim(ymaxswe);
        end

% RMSE vs n
    for cc = 1:length(n);
    for d = 1:8 
       den = DenOpt{d};
       stack_rmsen.(glacier)(d,cc) = subsetRmseLR(cc).(type).(den).(glacier);
    end
    end
subplot(3,3,6)
    if c==1
    plot(repmat(n,8,1),stack_rmsen.(glacier)(:,1:c),'.'); hold on
    else
    plot(n,stack_rmsen.(glacier)(:,1:c)','LineWidth',2); hold on
    end
    plot([0, max(Nfull)],[mean(rmseLR(g,:)), mean(rmseLR(g,:))],'k--')
        xlabel('Sample size'); ylabel('RMSE (m w.e.)')
        xlim([0, max(Nfull)])
        if c == 9; 
        ymaxrmse = ylim;
        else
        ylim(ymaxrmse);
        end

% Elevation Coeff vs n
 
    for cc = 1:length(n);
    for d = 1:8 
       den = DenOpt{d};
       stackC.(glacier)(d,cc) = (subsetLR(cc).(type).(den).coeff{1,g});
    end
    end

subplot(3,3,9)
    if c==1
    plot(repmat(n,8,1),stackC.(glacier)(:,1:c),'.'); hold on
    else
    plot(n,stackC.(glacier)(:,1:c)','LineWidth',2); hold on
    end
    plot([0, max(Nfull)],[mean(boxMLR.(glacier){1,:}), mean(boxMLR.(glacier){1,:})],'k--')
        xlabel('Sample size'); ylabel('\beta_z')
        xlim([0, max(Nfull)])
        if c == 9; 
        ymaxcoeff = ylim;
        else
        ylim(ymaxcoeff);
        end

%
     fig.PaperUnits = 'pixels'; fig.PaperPosition = [0 0 570 450];
 print(['/home/glaciology1/Documents/Data/Plots/SubsetAnimation/', glacier,type,num2str(c)],'-dpng','-r0');

end

%% Sitch together
close all
% load the images
 cd('/home/glaciology1/Documents/Data/Plots/SubsetAnimation/');
 fname = dir([glacier,type, '*.png']); images = cell(size(fname));
% fname = [fname(1); fname(3:end); fname(2)];
 for i = 1:length(fname);
 images{i} = imread( fname(i).name );
 end
 
 % create the video writer with 1 fps
 writerObj = VideoWriter(['Subset_',type,glacier,'.avi']);
 writerObj.FrameRate = 1;

 % set the seconds per image
 secsPerImage = [5 10 15];

 % open the video writer
 open(writerObj);

 % write the frames to the video
 for u=1:length(images)
     % convert the image to a frame
     frame = im2frame(images{u});

     %for v=1:secsPerImage(u) 
         writeVideo(writerObj, frame);
     %end
 end

 % close the writer object
 close(writerObj);
 
end