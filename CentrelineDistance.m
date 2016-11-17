%% Import data from text file.
centreline = csvread('centrelinepoints.csv');

%% Calculate distance from centre line

%Sampled points
glacier = {'G4','G2','G13'};
for i = 1:3
        G = char(glacier(i));

    Ecentre = repmat(centreline(:,1)',length(SWE(i).utm(:,1)),1); %repeating pit easting
    Ncentre = repmat(centreline(:,2)',length(SWE(i).utm(:,1)),1); %repeating pit northing
    Eloc = repmat(SWE(i).utm(:,1),1,length(centreline(:,1))); %repeating locations easting
    Nloc = repmat(SWE(i).utm(:,2),1,length(centreline(:,1))); %repeating locations northing

    X = Eloc-Ecentre; Y = Nloc-Ncentre; %separation between locations and pits
    distance = sqrt(X.^2+Y.^2);
    topo_sampled.(G).centreD = min(distance,[],2);
end
%scatter(SWE(i).utm(:,1), SWE(i).utm(:,2), 10, topo_sampled.G4.centreD,'filled'); 


%% Raster distance from centreline

%Create matrix with x and y locations of raster cells
corner = [593915 6739280; 598000 6751320; 602040 6760040]; %NW corner of raster

for i = 1:3
    G = char(glacier(i));    
    a = topo_full.(G).Sx;    a = flipud(a);     rasterSize = size(a);
    
    rasterX = corner(i,1):40:(40*rasterSize(1,2)+corner(i,1)-1);
    rasterX = repmat(rasterX, rasterSize(1,1),1);    
    rasterX(isnan(a)) = NaN;    rasterX = rasterX(:);

    rasterY = corner(i,2):40:(40*rasterSize(1,1)+corner(i,2)-1);
    rasterY = repmat(rasterY', 1, rasterSize(1,2));  
    rasterY(isnan(a)) = NaN;    rasterY = rasterY(:);

%     plot(rasterX,rasterY,'.'); hold on
%     plot(centreline(:,1),centreline(:,2),'k.');

Ecentre = repmat(centreline(:,1)',length(rasterX),1); %repeating centreline easting
Ncentre = repmat(centreline(:,2)',length(rasterY),1); %repeating centreline northing
Eloc = repmat(rasterX,1,length(centreline(:,1))); %repeating locations easting
Nloc = repmat(rasterY,1,length(centreline(:,1))); %repeating locations northing
    
X = Eloc-Ecentre; Y = Nloc-Ncentre; %separation between locations and pits
distance = sqrt(X.^2+Y.^2);
topo_full.(G).centreD = flipud(reshape(nanmin(distance,[],2),rasterSize));

clear E* N* raster* a
end