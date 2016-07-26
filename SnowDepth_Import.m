%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing Snow Depth Data
%       This script imports snowdepth data from the field data file for
%       transects, zigzags, and SWE. It also categorizes the group data

%       Alexandra Pulwicki July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing data

clear all
run MeasurementLocations.m  %This program determines the easting and northing of measurements

%Import data from excel sheet (change path based on computer)
    %[SD1, SD1text, SD1raw] = xlsread('/Volumes/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SD#1','A1:O833'); %Import values, text, and raw data from first sheet
    %[SD1, SD1text, SD1raw] = xlsread('/media/glaciology1/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SD#1','A1:O833'); %Import values, text, and raw data from first sheet
    [SD1, SD1text, SD1raw] = xlsread('FieldDataRevisedAP.xlsx','SD#1','A1:O833'); %Import values, text, and raw data from first sheet
    SD1text(1,:) = []; %Remaining text in sheet 1
    SD1(:,10) = []; %Remove comments column (not sure why it imports here)
    %[SD2, SD2text, SD2raw] = xlsread('/Volumes/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SD#2','A1:O832');
    %[SD2, SD2text, SD2raw] = xlsread('/media/glaciology1/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SD#2','A1:O832');
    [SD2, SD2text, SD2raw] = xlsread('FieldDataRevisedAP.xlsx','SD#2','A1:O832');
    SD2text(1,:) = [];
    %[SD3, SD3text, SD3raw] = xlsread('/Volumes/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SD#3','A1:O675'); 
    %[SD3, SD3text, SD3raw] = xlsread('/media/glaciology1/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SD#3','A1:O675'); 
    [SD3, SD3text, SD3raw] = xlsread('FieldDataRevisedAP.xlsx','SD#3','A1:O675'); 
    SD3text(1,:) = [];
    %[ZZ, ZZtext, ZZraw] = xlsread('/Volumes/GlacierAlex/Data/FieldDataRevisedAP.xlsx','ZigZag','A1:H1653'); 
    %[ZZ, ZZtext, ZZraw] = xlsread('/media/glaciology1/GlacierAlex/Data/FieldDataRevisedAP.xlsx','ZigZag','A1:H1653'); 
    [ZZ, ZZtext, ZZraw] = xlsread('FieldDataRevisedAP.xlsx','ZigZag','A1:H1653'); 
    ZZtext(1,:) = [];
    %[ExtraSD, ExtraSDtext, ExtraSDraw] = xlsread('/Volumes/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SWEDepth','A1:AU38'); 
    %[ExtraSD, ExtraSDtext, ExtraSDraw] = xlsread('/media/glaciology1/GlacierAlex/Data/FieldDataRevisedAP.xlsx','SWEDepth','A1:AU38'); 
    [ExtraSD, ExtraSDtext, ExtraSDraw] = xlsread('FieldDataRevisedAP.xlsx','SWEDepth','A1:AU38'); 
    ExtraSDtext(1,:) = [];

%Import vertex waypoints
    %[~,~,Vertex_cord] = xlsread('/Volumes/GlacierAlex/QGIS/Donjek Glaciers/Sampling/zigzag_corners_utm.xls');
    %[~,~,Vertex_cord] = xlsread('/media/glaciology1/GlacierAlex/QGIS/Donjek Glaciers/Sampling/zigzag_corners_utm.xls');
    [~,~,Vertex_cord] = xlsread('/home/glaciology1/Documents/QGIS/Donjek Glaciers/Sampling/zigzag_corners_utm.xls');

%% Categorizing data

%Categorize quality and situational data
    %Converting the text to categories means that it can be used to group
    %data for plotting and stats analysis in the future
%Transect data
    SD1_Q = [categorical(SD1(:,3)),categorical(SD1(:,5)),...    %Quality of data (1 for good, 0 for question mark)
        categorical(SD1(:,7)),categorical(SD1(:,9))];           
    SD1_Glacier = categorical(SD1text(:,12));                   %Glacier (G13, G02, G04)
    SD1_Person = categorical(SD1text(:,13));                    %Person that took the measurement (AP, CA, GF, AC)
    SD1_Pattern = categorical(SD1text(:,14));                   %Transect patterns that the point is a part of
    SD1_Date = categorical(SD1text(:,15));                      %Date meaurement was taken
    SD1_Book = categorical(SD1text(:,11));                      %Book/order for measurement
SD1 = [SD1(:,2),SD1(:,4),SD1(:,6),SD1(:,8),SD1(:,1)];       %Depth 1, 2, 3, 4, WP#

    SD2_Q = [categorical(SD2(:,3)),categorical(SD2(:,5)),...
        categorical(SD2(:,7)),categorical(SD2(:,9))];
    SD2_Glacier = categorical(SD2text(:,12));
    SD2_Person = categorical(SD2text(:,13));
    SD2_Pattern = categorical(SD2text(:,14));
    SD2_Date = categorical(SD2text(:,15));
    SD2_Book = categorical(SD2text(:,11));
SD2 = [SD2(:,2),SD2(:,4),SD2(:,6),SD2(:,8),SD2(:,1)];       

    SD3_Q = [categorical(SD3(:,3)),categorical(SD3(:,5)),...
        categorical(SD3(:,7)),categorical(SD3(:,9))];
    SD3_Glacier = categorical(SD3text(:,12));
    SD3_Person = categorical(SD3text(:,13));
    SD3_Pattern = categorical(SD3text(:,14));
    SD3_Date = categorical(SD3text(:,15));
    SD3_Book = categorical(SD3text(:,11));
SD3 = [SD3(:,2),SD3(:,4),SD3(:,6),SD3(:,8),SD3(:,1)];       
%Zigzag data
    ZZ_Glacier = categorical(ZZtext(:,1));                      %Glacier (G13, G02, G04)
    ZZ_Zone = categorical(ZZtext(:,2));                         %Zone label
    ZZ_Vertex = categorical(ZZtext(:,3));                       %Reference vertex
    ZZ_Q = categorical(ZZ(:,3));                                %Quality of data (1 for good, 0 for question mark)
    ZZ_Person = categorical(ZZtext(:,8));                       %Person that took the measurement (AP, CA, GF, AC)
    ZZ_Date = categorical(ZZtext(:,7));                         %Date meaurement was taken
    C(:) = {'ZZ'};
    ZZ_Book = categorical(C);
ZZ = [zeros(length(ZZ),1),ZZ];                              %[zeros(will be WP#), distance data, depth, quality]
%Depth data from SWE and snowpits
    ExtraSD_Glacier = categorical(ExtraSDtext(:,44));           %%Glacier (G13, G02, G04)
    C = cell(size(ExtraSD_Glacier));
    C(:) = {'Extra'};
    ExtraSD_Pattern = categorical(C);
    ExtraSD_Person = categorical(C);
    ExtraSD_Book = categorical(C);
    C(:) = {'1'};
    ExtraSD_Q = categorical(C);
        clear C
               
%% Working with transect data

% Getting easting and northing for each data point (columns 6 and 7) from
% the vector 'closest' from the MeasurementLocation.m script
for i = 1:size(SD1,1) %for all the waypoints in the SD1 matrix
    ind = find(SD1(i,5)==floor(closest(:,3))); % get index of corresponding waypoint in 'closest' (there will be three because 4.1, 4.2, 4.3)
    SD1(i,6:7) = closest(ind(1),1:2); %assign the first easting and northing from 'closest' to column 6 and 7 of SD
end

for i = 1:size(SD2,1)
    ind = find(SD2(i,5)==floor(closest(:,3)));
    SD2(i,6:7) = closest(ind(2),1:2);
end

for i = 1:size(SD3,1)
    ind = find(SD3(i,5)==floor(closest(:,3)));
    SD3(i,6:7) = closest(ind(3),1:2);
end

clear ind i closest

% Getting all waypoints in SD1, SD2, and SD3 data (empty for ones not 
% measured). This allows SD matrices to be the same size and have aligned
% WPs and indices
for i=1:max(SD1(:,5))-3 %use WP range so that gaps can be filled
   if i+3~=SD1(i,5) %determine if WP is missing from the data (first WP is 4 so each WP# = index+3)
       SD1 = [SD1(1:i-1,:); nan,nan,nan,nan,i+3,nan,nan ;SD1(i:end,:)]; %insert nan row with WP# in SD1
       SD1_Date = [SD1_Date(1:i-1,1); '<undefined>' ;SD1_Date(i:end,1)]; %insert and undefinied into categorical arrays
       SD1_Glacier = [SD1_Glacier(1:i-1,1); '<undefined>' ;SD1_Glacier(i:end,1)];
       SD1_Pattern = [SD1_Pattern(1:i-1,1); '<undefined>' ;SD1_Pattern(i:end,1)];
       SD1_Person = [SD1_Person(1:i-1,1); '<undefined>' ;SD1_Person(i:end,1)];
       SD1_Book = [SD1_Book(1:i-1,1); '<undefined>' ;SD1_Book(i:end,1)];
       SD1_Q = [SD1_Q(1:i-1,1:4);{'<undefined>', '<undefined>','<undefined>','<undefined>'};SD1_Q(i:end,1:4)];
       SD1raw = [SD1raw(1:i-1,:); repmat({NaN},1,15) ;SD1raw(i:end,:)]; %insert nan rows into cell arrays
       SD1text = [SD1text(1:i-1,:); repmat({'None'},1,15) ;SD1text(i:end,:)];
   end
   if i+3~=SD2(i,5)
       SD2 = [SD2(1:i-1,:); nan,nan,nan,nan,i+3,nan,nan ;SD2(i:end,:)];
       SD2_Date = [SD2_Date(1:i-1,1); '<undefined>' ;SD2_Date(i:end,1)];
       SD2_Glacier = [SD2_Glacier(1:i-1,1); '<undefined>' ;SD2_Glacier(i:end,1)];
       SD2_Pattern = [SD2_Pattern(1:i-1,1); '<undefined>' ;SD2_Pattern(i:end,1)];
       SD2_Person = [SD2_Person(1:i-1,1); '<undefined>' ;SD2_Person(i:end,1)];
       SD2_Book = [SD2_Book(1:i-1,1); '<undefined>' ;SD2_Book(i:end,1)];
       SD2_Q = [SD2_Q(1:i-1,1:4);{'<undefined>', '<undefined>','<undefined>','<undefined>'};SD2_Q(i:end,1:4)];
       SD2raw = [SD2raw(1:i-1,:); repmat({NaN},1,15) ;SD2raw(i:end,:)];
       SD2text = [SD2text(1:i-1,:); repmat({'None'},1,15) ;SD2text(i:end,:)];
   end
   if i+3~=SD3(i,5)
       SD3 = [SD3(1:i-1,:); nan,nan,nan,nan,i+3,nan,nan ;SD3(i:end,:)];
       SD3_Date = [SD3_Date(1:i-1,1); '<undefined>' ;SD3_Date(i:end,1)];
       SD3_Glacier = [SD3_Glacier(1:i-1,1); '<undefined>' ;SD3_Glacier(i:end,1)];
       SD3_Pattern = [SD3_Pattern(1:i-1,1); '<undefined>' ;SD3_Pattern(i:end,1)];
       SD3_Person = [SD3_Person(1:i-1,1); '<undefined>' ;SD3_Person(i:end,1)];
       SD3_Book = [SD3_Book(1:i-1,1); '<undefined>' ;SD3_Book(i:end,1)];
       SD3_Q = [SD3_Q(1:i-1,1:4);{'<undefined>', '<undefined>','<undefined>','<undefined>'};SD3_Q(i:end,1:4)];
       SD3raw = [SD3raw(1:i-1,:); repmat({NaN},1,15) ;SD3raw(i:end,:)];
       SD3text = [SD3text(1:i-1,:); repmat({'None'},1,15) ;SD3text(i:end,:)];
   end
   
end
clear i 

SD1text(cellfun(@isempty,SD1text(:,10)),10)={'None'};
SD2text(cellfun(@isempty,SD2text(:,10)),10)={'None'};
SD3text(cellfun(@isempty,SD3text(:,10)),10)={'None'};
ExtraSDtext(cellfun(@isempty,ExtraSDtext(:,43)),43)={'None'};


field1 = 'raw';         value1 = {SD1raw, SD2raw, SD3raw, ExtraSDraw};
field2 = 'depth';       value2 = {SD1, SD2, SD3, ExtraSD};
field3 = 'glacier';     value3 = {SD1_Glacier, SD2_Glacier, SD3_Glacier, ExtraSD_Glacier};
field4 = 'pattern';     value4 = {SD1_Pattern, SD2_Pattern, SD3_Pattern, ExtraSD_Pattern};
field5 = 'person';      value5 = {SD1_Person, SD2_Person, SD3_Person, ExtraSD_Person};
field6 = 'Q';           value6 = {SD1_Q, SD2_Q, SD3_Q, ExtraSD_Q};
field7 = 'book';        value7 = {SD1_Book, SD2_Book, SD3_Book, ExtraSD_Book};
field8 = 'comments';    value8 = {SD1text(:,10), SD2text(:,10), SD3text(:,10), ExtraSDtext(:,43)};

SD = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
    field5,value5,field6,value6,field7,value7,field8,value8);
clear value* field* SD1* SD2* SD3* ExtraSD*


%% Working with zigzag data

field1 = 'raw';         value1 = {ZZraw};
field2 = 'depth';       value2 = {ZZ};
field3 = 'vertex';      value3 = {ZZ_Vertex};
field4 = 'zone';        value4 = {ZZ_Zone};
field5 = 'glacier';     value5 = {ZZ_Glacier};
field6 = 'person';      value6 = {ZZ_Person};
field7 = 'Q';           value7 = {ZZ_Q};
field8 = 'book';        value8 = {ZZ_Book};

ZZ = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
    field5,value5,field6,value6,field7,value7,field8,value8);
clear value* field* ZZ_* ZZraw ZZtext
