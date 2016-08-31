%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Calculating the locations of snow depth measurements
%Alexandra Pulwicki
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Import data from excel sheet
    [GPS, ~, GPSraw] = xlsread('GlacierWP_UTM.xlsx','A1:G845'); %Import values, text, and raw data from first sheet
    GPS(:,6) = [];

% Allocate imported array to column variable names
    easting = GPS(:, 1);
    northing = GPS(:, 2);
    lat = GPS(:, 3);
    long = GPS(:, 4);
    gps_elev = GPS(:, 5);
    time = GPSraw(1:end,6);
    name = GPS(:, 6);
    
%% Creating imaginary waypoints for the inital waypoint of each shape

% Waypoint number of the start of each shape
    starters = [4,21,72,127,159,185,208,...
       223,276,314,344,355,371,519,...
       529,571,660,678,714,745,812];

for i = 1:length(starters) %for each start waypoint
    index = find(name==starters(i));

    [~,m,b] = regression(easting(index:index+1)',northing(index:index+1)'); %Generates stright line between start WP and subsequent WP
    x = easting(index)+40; %Easting of imaginary WP
    y = m*x+b; %Northing of imaginary WP

        %This part changes the location of the imaginary WP so that is it
        %behind the start WP (not between the start and subsequent WPs)
    if  EuclideanDistance(easting(index,1),northing(index,1),x,y) > EuclideanDistance(easting(index+1,1),northing(index+1,1),x,y) %The distance between the imaginary and start WP should be smaller that that between the imaginary and subsequent WP
        x = easting(index)-40; %New easting of imaginary WP
        y = m*x+b; %New northing of imaginary WP
    end
    
        %Insert the imaginary WP into the data
    easting = [easting(1:index-1);x;easting(index:end)];
    northing = [northing(1:index-1);y;northing(index:end)];
    name = [name(1:index-1);str2double(strcat(num2str(name(index)),'000'));name(index:end)]; %The name of the imaginary WP is 'OldWP#' and then 000

end

%% Location of actual depth measurements based on GPS coordinates
    closest = zeros(3,1); %Create matrix for the points between WPs (aka actual measurement locations)

for j = 1:length(easting)-1 %for each WP

    if name(j+1)<1000 %ensures that two subsequent points are from the same shape

        % Create a series of points between two GPS locations 
        line = linspaceNDim([easting(j,1),northing(j,1)],[easting(j+1,1),northing(j+1,1)],1000);

        % Find the distance between created points and the original GPS location
        distinterp = zeros(length(line),1);
        for i = 1:length(line)
            distinterp(i,1) = EuclideanDistance(line(1,i),line(2,i),easting(j+1,1),northing(j+1,1));
        end

        % Selects points closest to the chosen distance from the GPS
        % The first part is for when there were only two probers (G02 LH &
        % LC). The second is for when there were three probers (remainder).
        if name(j)>370 && name(j)<529 || name(j) == 3710 || name(j) == 5190 %WPs when only two probers
            val = [10 20]; %Chosen distances from GPS of probers
            [~, index] = min(abs(distinterp-val(1,1))); %Finds the index of the point closest to the first chosen distances
            closest(:,end+1) = [line(:,index);(name(j+1)+0.1)]; %Enters the easting, northing, and GPS WP into variable 'closest' for SD1 book
            [~, index] = min(abs(distinterp-val(1,2)));
            closest(:,end+1) = [line(:,index);(name(j+1)+0.2)]; %Enters the easting, northing, and GPS WP into variable 'closest' for SD2 book
        else %WPs when three probers
            val = [10 20 30]; %Chosen distance from GPS of probers
            [~, index] = min(abs(distinterp-val(1,1)));
            closest(:,end+1) = [line(:,index);(name(j+1)+0.1)]; %SD1 (eg. WP = 4.1)
            [~, index] = min(abs(distinterp-val(1,2)));
            closest(:,end+1) = [line(:,index);(name(j+1)+0.2)]; %SD2 (eg. WP = 4.2)
            [~, index] = min(abs(distinterp-val(1,3)));
            closest(:,end+1) = [line(:,index);(name(j+1)+0.3)]; %SD3 (eg. WP = 4.3)
        end
    end
    
end
    clear distinterp index i j line x y m b 

    clear easting northing lat long name starters time val GPS GPSraw

closest = closest(:,2:end)'; %Transpose matrix and remove initial 0

%dlmwrite('/home/glaciology1/Documents/QGIS/Data/TestPoints.csv',closest, 'delimiter', ',', 'precision', 9); %Write matrix with new waypoints to csv file for QGIS

