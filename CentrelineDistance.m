% Calculating the distance of all raster points from the centreline 
%       This script finds the smallest distance between all sampled point 
%       and all points in a raster from the glacier centreline (drawn by 
%       hand in QGIS). 
%
%       Inputs:         centrelinepoints.csv (from QGIS)
%       Outputs:        topo_full.(G).centreD

%       Alexandra Pulwicki  Created: October 2016
%                           Updated: November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Import centreline points
centreline = csvread('centrelinepoints.csv');

%% Distance of sampled points from centre line

for i = 1:3
    G = char(options.glacier(i));

    Ecentre = repmat(centreline(:,1)',length(SWE(i).utm(:,1)),1); %repeating centreline easting
    Ncentre = repmat(centreline(:,2)',length(SWE(i).utm(:,1)),1); %repeating centreline northing
    Eloc = repmat(SWE(i).utm(:,1),1,length(centreline(:,1))); %repeating sampled locations easting
    Nloc = repmat(SWE(i).utm(:,2),1,length(centreline(:,1))); %repeating sampled locations northing

    X = Eloc-Ecentre; Y = Nloc-Ncentre; %separation between locations and centreline
    distance = sqrt(X.^2+Y.^2); %calculate distance
    topo_sampled.(G).centreD = min(distance,[],2); %chose the minimum distance
end
%scatter(SWE(i).utm(:,1), SWE(i).utm(:,2), 10, topo_sampled.G4.centreD,'filled'); 

%% Raster distance from centreline

%Create matrix with x and y locations of raster cells
corner = [593915 6739280; 598000 6751320; 602040 6760040]; %Coordinates of NW corner of rasters

for i = 1:3
    G = char(options.glacier(i));    
    a = topo_full.(G).Sx;    a = flipud(a);     rasterSize = size(a); %get example raster for NaN locations and it's size
    
    rasterX = corner(i,1):40:(40*rasterSize(1,2)+corner(i,1)-1); %create new raster with easting locatiosn of cells
    rasterX = repmat(rasterX, rasterSize(1,1),1);    
    rasterX(isnan(a)) = NaN;    rasterX = rasterX(:); %put NaNs where they need to be

    rasterY = corner(i,2):40:(40*rasterSize(1,1)+corner(i,2)-1); %create new raster with northing locatiosn of cells
    rasterY = repmat(rasterY', 1, rasterSize(1,2));  
    rasterY(isnan(a)) = NaN;    rasterY = rasterY(:); %put NaNs where they need to be

%     plot(rasterX,rasterY,'.'); hold on
%     plot(centreline(:,1),centreline(:,2),'k.');

Ecentre = repmat(centreline(:,1)',length(rasterX),1); %repeating centreline easting
Ncentre = repmat(centreline(:,2)',length(rasterY),1); %repeating centreline northing
Eloc = repmat(rasterX,1,length(centreline(:,1))); %repeating raster cell locations easting
Nloc = repmat(rasterY,1,length(centreline(:,1))); %repeating raster cell locations northing
    
X = Eloc-Ecentre; Y = Nloc-Ncentre; %separation between locations and cells
distance = sqrt(X.^2+Y.^2); %calculate distance
topo_full.(G).centreD = flipud(reshape(nanmin(distance,[],2),rasterSize)); %chose the minimum distance

clear E* N* raster* a
end