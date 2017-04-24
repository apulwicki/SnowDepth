function [ dataout ] = KrigingR_G( data )
%UNTITLED2 Summary of this function goes here
global options

for g = 1:3
    glacier = options.glacier{g};

%   Detailed explanation goes here
% Dice Kriging -> residuals
    clear utm
    utm(:,1)    = data.(glacier)(:,2)-min(options.rig.(glacier)(:,1));
    utm(:,2)    = data.(glacier)(:,3)-min(options.rig.(glacier)(:,2));
    res         = data.(glacier)(:,1);
    sizexy      = options.mapsize(g,:);

cd Kriging

    save('residuals.mat','res','utm','sizexy')

    %Run Dice Kriging in R
    !R CMD BATCH DiceKriging.R
    %!/usr/local/bin/R CMD BATCH DiceKriging.R

    % load kriged data
    load kriging.mat
    
     %Model params
    dataout.Model   =  model;
%     dataout.LOO     =  LOO.mean;

     %Assign to structure
    dataout.(glacier)    = flipud(pred);
%     dataout.lower95 = flipud(lower95);
%     dataout.upper95 = flipud(upper95);
        
     %Set glacier boundries
    dataout.(glacier)(options.mapNaN.(glacier))    = NaN;
%     dataout.lower95(options.mapNaN.(glacier)) = NaN;
%     dataout.upper95(options.mapNaN.(glacier)) = NaN;
        
cd ..

end
end

