%Import data
    filename = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/SianCampWeather.txt';
    delimiter = '\t';
    formatSpec = '%q%*q%*q%*q%*q%*q%f%*s%[^\n\r]';
    fileID = fopen(filename,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'EmptyValue' ,NaN, 'ReturnOnError', false);
    fclose(fileID);
    dataArray(2) = cellfun(@(x) num2cell(x), dataArray(2), 'UniformOutput', false);
SianCampWeather = [dataArray{1:end-1}];
    clearvars filename delimiter formatSpec fileID dataArray ans;
    
   Weather = zeros(length(SianCampWeather),5);
for i = 1:length(SianCampWeather)
   T = SianCampWeather{i,1}; 
   Weather(i,1) = str2double(T(1,4:5)); %day
   Weather(i,2) = str2double(T(1,1:2)); %month
   Weather(i,3) = str2double(T(1,7:8)); %year
   Weather(i,4) = str2double(T(1,10:11)); %hour
   Weather(i,5) = cell2mat(SianCampWeather(i,2)); %temp
end

%Scale to G13 height
    G13Weather = Weather;
    lapseR = -6.5; %K/km
    stationElev = 3050/1000;
    G13Elev = 2054/1000;
G13Weather(:,5) = Weather(:,5)+ lapseR*(G13Elev-stationElev);
%% 

    PosDay = Weather(Weather(:,5)>0,:);
    DDfact = 4; %mm/day/K
    
Melt = DDfact*1/24*PosDay(:,5); 
TotalMelt = sum(Melt)/1000; %melt in m w.e.