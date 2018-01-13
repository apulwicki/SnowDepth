
   %Loads input data: point-scale b_w, topographic parameters from sampled
    %cells, and options
    load TopoSWE.mat fullSWE topo_sampled topo_full
    run OPTIONS

    %Generates gridcell-averaged b_w values for input to the OK algorithm
    for d = 1:8;    den = options.DenOpt{d};
    [ inputSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled);
    end


for g = 1:3
    cd Kriging
    
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
    
    cd ..
    
    predELEV(predELEV==predELEV(1,1)) = NaN;
    UKelev.(glacier) = predELEV;
    
    subplot(1,3,g)
        imagesc(predELEV); colorbar
        
end