%% Import data from text file.
filename = '/home/glaciology1/Documents/Data/SnowDepth/centrelinepoints.csv';
delimiter = ',';
formatSpec = '%s%s%[^\n\r]';
fileID = fopen(filename,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,  'ReturnOnError', false);
fclose(fileID);
raw = repmat({''},length(dataArray{1}),length(dataArray)-1);
for col=1:length(dataArray)-1
    raw(1:length(dataArray{col}),col) = dataArray{col};
end
numericData = NaN(size(dataArray{1},1),size(dataArray,2));
for col=[1,2]
    % Converts strings in the input cell array to numbers. Replaced non-numeric
    % strings with NaN.
    rawData = dataArray{col};
    for row=1:size(rawData, 1);
        % Create a regular expression to detect and remove non-numeric prefixes and
        % suffixes.
        regexstr = '(?<prefix>.*?)(?<numbers>([-]*(\d+[\,]*)+[\.]{0,1}\d*[eEdD]{0,1}[-+]*\d*[i]{0,1})|([-]*(\d+[\,]*)*[\.]{1,1}\d+[eEdD]{0,1}[-+]*\d*[i]{0,1}))(?<suffix>.*)';
        try
            result = regexp(rawData{row}, regexstr, 'names');
            numbers = result.numbers;
            
            % Detected commas in non-thousand locations.
            invalidThousandsSeparator = false;
            if any(numbers==',');
                thousandsRegExp = '^\d+?(\,\d{3})*\.{0,1}\d*$';
                if isempty(regexp(numbers, thousandsRegExp, 'once'));
                    numbers = NaN;
                    invalidThousandsSeparator = true;
                end
            end
            % Convert numeric strings to numbers.
            if ~invalidThousandsSeparator;
                numbers = textscan(strrep(numbers, ',', ''), '%f');
                numericData(row, col) = numbers{1};
                raw{row, col} = numbers{1};
            end
        catch me
        end
    end
end
R = cellfun(@(x) ~isnumeric(x) && ~islogical(x),raw); % Find non-numeric cells
raw(R) = {NaN}; % Replace non-numeric cells

centreline = cell2mat(raw); centreline(1,:) = [];

clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;


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