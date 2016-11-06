filename = 'centrelinepoints.csv';
delimiter = ',';
formatSpec = '%q%q%[^\n\r]';
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
X1 = cell2mat(raw(:, 1));
Y1 = cell2mat(raw(:, 2));
clearvars filename delimiter formatSpec fileID dataArray ans raw col numericData rawData row regexstr result numbers invalidThousandsSeparator thousandsRegExp me R;

centreline = [X1 Y1]; centreline(1,:) = []; clear X1 Y1
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
    min_dist.(G) = min(distance,[],2);
end


%Raster points

%Create matrix with x and y locations of raster cells
corner = [centreline(1,1) centreline(1,2)]; %NW corner of raster
rasterSize = size(topo_full.G4.Sx);

rasterX = corner(1,1):40:(40*rasterSize(1,2)+corner(1,1)-1);
rasterX = repmat(rasterX, rasterSize(1,1),1);   rasterX = rasterX(:);
rasterY = [corner(1,2):40:(40*rasterSize(1,1)+corner(1,2)-1)]';
rasterY = repmat(rasterY, 1, rasterSize(1,2));  rasterY = rasterY(:);

Ecentre = repmat(centreline(:,1)',length(rasterX),1); %repeating centreline easting
Ncentre = repmat(centreline(:,2)',length(rasterY),1); %repeating centreline northing
Eloc = repmat(rasterX,1,length(centreline(:,1))); %repeating locations easting
Nloc = repmat(rasterY,1,length(centreline(:,1))); %repeating locations northing
    
X = Eloc-Ecentre; Y = Nloc-Ncentre; %separation between locations and pits
distance = sqrt(X.^2+Y.^2);
min_dist.(G) = min(distance,[],2);