function [ dataout ] = KrigingR( data, coord, glacier)
%UNTITLED2 Summary of this function goes here
global options

if strcmp(glacier, 'G4')
    g = 1;
elseif strcmp(glacier, 'G2')
    g = 2;
elseif strcmp(glacier, 'G13')
    g = 3;
end

%   Detailed explanation goes here
% Dice Kriging -> residuals

    utm(:,1)    = coord(:,1)-min(options.rig.(glacier)(:,1));
    utm(:,2)    = coord(:,2)-min(options.rig.(glacier)(:,2));
    res         = data;
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
    dataout.LOO     = LOO.mean;

     %Assign to structure
    dataout.pred    = flipud(pred);
    dataout.lower95 = flipud(lower95);
    dataout.upper95 = flipud(upper95);
        
     %Set glacier boundries
    dataout.pred(options.mapNaN.(glacier))    = NaN;
    dataout.lower95(options.mapNaN.(glacier)) = NaN;
    dataout.upper95(options.mapNaN.(glacier)) = NaN;
        
cd ..

end

