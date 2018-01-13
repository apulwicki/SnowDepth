for g = 1:3
    glacier = options.glacier{g};

%   Detailed explanation goes here
% Dice Kriging -> residuals
    clear utm
    data        = inputSWE.S2;
    
    utm         = data.(glacier)(:,4);
    res         = data.(glacier)(:,1);
    sizexy      = options.mapsize(g,:);
    EE          = topo_full.(glacier).elevation;
    EE(isnan(EE))=0;

    save('residuals.mat','res','utm','sizexy','EE')

    %Run Dice Kriging in R
    %!R CMD BATCH DiceKriging.R
    !/usr/local/bin/R CMD BATCH UKtesting.R

    % load kriged data
    load UKelevTest.mat
    
    predELEV(predELEV==predELEV(1,1)) = NaN;
    UKelev.(glacier) = predELEV;
    
    subplot(1,3,g)
        imagesc(predELEV); colorbar
        
end