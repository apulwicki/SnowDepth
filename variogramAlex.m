function [ d ] = variogramAlex(zz, lag, maxlag_in)
%Calculates the variogam values for input data
%   This function calculates the variance between two points and their
%   distance apart. It then bins the pairs of points into specified lags
%   and calculates the average variance (semi-variance). The results are
%   returned as a strctured array with semi-variance, binned distance, and
%   number of pairs. Results are plotted. 

%   Lag is speficied by the user (lag) and the lag tolerance is set as half 
%   the lag. For example, if the lag is set at 10 then the distances will be
%   binned between 0-10, 10-20 etc, and the points plotted at 5, 15 etc.
%   Maximum lag distance can be specified but can be set to 'default' where
%   it is determined to be half the domain distance. 

%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Check formating of input data - converts to [value x y]
    if isstruct(zz) %identifies transect data (fat format)
        data = [nanmean(zz(5).depth(:,1:4),2), zz(5).depth(:,6:7)]; %calculates mean depth
    else %identifes zigzag data
        data = zz;
    end
    
% Calculate distance and variance between pairs of points 
    z = data(:,1); x = data(:,2); y = data(:,3);
    X = meshgrid(x); Y = meshgrid(y); Z = meshgrid(z);
    
    distMtx = triu(sqrt((X-X').^2 + (Y-Y').^2));
    varMtx = triu((Z-Z').^2);
    
    % Set maxlag if not specified    
    if strcmp(maxlag_in,'default') %when 'default' is chosen
        maxlag = ceil(max(max(distMtx))/lag/2.5)*lag; %maxlag is half the maximum domain distance (rounded up to have integer number of lags)
    else
        maxlag = ceil(maxlag_in/lag)*lag;
    end
    
    bins = 0:lag:maxlag;    
    binCentres = bins(1:end-1)+(lag/2);
    distBin = zeros(length(bins)-1,1); semiVar = distBin; numEls = distBin;
    for i = 1:length(bins)-1
        ll = bins(i); ul = bins(i+1);
        logiEl = distMtx > ll & distMtx <= ul;
        
        numEls(i) = sum(sum(logiEl));
        
        distBin(i) = mean(mean(distMtx(logiEl)));
        
        varBin = sum(sum(varMtx(logiEl)));
        semiVar(i) = varBin/(2*numEls(i));
    end
       
% Create structure with output data
    d = struct('meanDist', distBin, 'val', semiVar, 'num', numEls, 'binCentre', binCentres');

end

