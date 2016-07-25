%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Importing Snow Depth Data
%       This script imports snowdepth data from the field data file for
%       transects, zigzags, and SWE. It also categorizes the group data

%       Alexandra Pulwicki July 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing data

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
ZZ = [zeros(length(ZZ),1),ZZ];                              %[zeros(will be WP#), distance data, depth, quality]
%Depth data from SWE and snowpits
    ExtraSD_Glacier = categorical(ExtraSDtext(:,44));           %%Glacier (G13, G02, G04)
    C = cell(size(ExtraSD_Glacier));
    C(:) = {'ZZ'};
    ExtraSD_Pattern = categorical(C);
    ExtraSD_Person = categorical(C);
    ExtraSD_Book = categorical(C);
    C(:) = {'1'};
    ExtraSD_Q = categorical(C);
        clear C
