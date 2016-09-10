% Importing Transect Snow Depth Data
%       This script imports snowdepth data from the field data file for
%       transects, zigzags, and SWE. It also categorizes the descriptive
%       part of the data. Then it creates a structured array for transect
%       data
%
%       Inputs:         Field Data ('FieldDataRevisedAP.xlsx')
%       Outputs:        Snowdepth structure (SD)
%                       Elevations from GPS WPs (gps_elev)

%       Alexandra Pulwicki  Created: July 2016
%                           Updated: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing data

%Import data from excel sheet (change path based on computer)
    [SD1, SD1text, SD1raw] = xlsread('FieldDataRevisedAP.xlsx','SD#1','A1:O833'); %Import values, text, and raw data from first sheet
    SD1text(1,:) = []; %Remaining text in sheet 1

    [SD2, SD2text, SD2raw] = xlsread('FieldDataRevisedAP.xlsx','SD#2','A1:O832');
    SD2text(1,:) = [];

    [SD3, SD3text, SD3raw] = xlsread('FieldDataRevisedAP.xlsx','SD#3','A1:O675'); 
    SD3text(1,:) = [];

    [ExtraSD, ExtraSDtext, ExtraSDraw] = xlsread('FieldDataRevisedAP.xlsx','SWEDepth','A1:AU38'); 
    ExtraSDtext(1,:) = [];

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

%Depth data from SWE and snowpits (aka ExtraSD)
    ExtraSD_Glacier = categorical(ExtraSDtext(:,44));           %Glacier (G13, G02, G04)
        WP_index = [9,20:28,34];
        for i = 1:length(WP_index)
            ExtraSDraw{WP_index(i),1}=num2str(ExtraSDraw{WP_index(i),1});
        end
    ExtraSD_Person = categorical(ExtraSDraw(2:end,1));           %Person (waypoint label)
        C = cell(size(ExtraSD_Glacier));                            
        C(:) = {'Extra'};
    ExtraSD_Pattern = categorical(C);                            %Pattern (all Extra)
    ExtraSD_Book = categorical(C);                               %Book (all Extra)
        C(:) = {'1'};
    ExtraSD_Q = categorical(C);                                  %Quality (all 1 = good)
        clear C
               
%% Working with transect data

% Getting easting and northing for each data point (columns 6 and 7) from
% the vector 'closest' from the MeasurementLocation.m script
    for i = 1:size(SD1,1) %for all the waypoints in the SD1 matrix
        ind = find(SD1(i,5)==floor(closest(:,3))); %get index of corresponding waypoint in 'closest' (there will be three because 4.1, 4.2, 4.3)
        SD1(i,6:7) = closest(ind(1),1:2); %assign the first easting and northing from 'closest' to column 6 and 7 of SD
    end

    for i = 1:size(SD2,1)
        ind = find(SD2(i,5)==floor(closest(:,3)));
        SD2(i,6:7) = closest(ind(2),1:2); %assign the second easting and northing from 'closest' to column 6 and 7 of SD
    end

    for i = 1:size(SD3,1)
        ind = find(SD3(i,5)==floor(closest(:,3)));
        SD3(i,6:7) = closest(ind(3),1:2); %assign the third easting and northing from 'closest' to column 6 and 7 of SD
    end
        clear ind i closest

% Getting all waypoints in SD1, SD2, and SD3 data (empty for ones not 
% measured). This allows SD matrices to be the same size and have aligned
% WPs and indices
for i=1:max(SD1(:,5))-3 %use WP range so that gaps can be filled
   if i+3~=SD1(i,5) %determine if WP is missing from the data (first WP is 4 so each WP# = index+3)
       SD1 = [SD1(1:i-1,:); nan,nan,nan,nan,i+3,nan,nan ;SD1(i:end,:)]; %insert nan row with WP# in SD1
       SD1_Date = [SD1_Date(1:i-1,1); '<undefined>' ;SD1_Date(i:end,1)]; %insert an undefinied into categorical arrays
       SD1_Glacier = [SD1_Glacier(1:i-1,1); '<undefined>' ;SD1_Glacier(i:end,1)];
       SD1_Pattern = [SD1_Pattern(1:i-1,1); '<undefined>' ;SD1_Pattern(i:end,1)];
       SD1_Person = [SD1_Person(1:i-1,1); '<undefined>' ;SD1_Person(i:end,1)];
       SD1_Book = [SD1_Book(1:i-1,1); '<undefined>' ;SD1_Book(i:end,1)];
       SD1_Q = [SD1_Q(1:i-1,1:4);{'<undefined>', '<undefined>','<undefined>','<undefined>'};SD1_Q(i:end,1:4)];
       SD1raw = [SD1raw(1:i-1,:); repmat({NaN},1,15) ;SD1raw(i:end,:)]; %insert nan rows into cell arrays
       SD1text = [SD1text(1:i-1,:); repmat({'None'},1,15) ;SD1text(i:end,:)]; %insert 'None' text into text array
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

    SD1text(cellfun(@isempty,SD1text(:,10)),10)={'None'};       %fills in the empty cells in the comments columns with 'None'
    SD1comments = categorical(SD1text(:,10));                   %Categorizes comments written

    SD2text(cellfun(@isempty,SD2text(:,10)),10)={'None'};
    SD2comments = categorical(SD2text(:,10));                   

    SD3text(cellfun(@isempty,SD3text(:,10)),10)={'None'};
    SD3comments = categorical(SD3text(:,10));                  
    
    ExtraSDtext(cellfun(@isempty,ExtraSDtext(:,43)),43)={'None'};
    ExtraSDcomments = categorical(ExtraSDtext(:,43));                  

%% Creating the structure for transect snowdepth data -> SD
    %Field values are names for columns, values are the cell arrays within the
    %structure. Original data arrangement is retained. SD1 is in row 1, 
    %SD2 in row 2, SD3 in row 3, ExtraSD in row 4. Example of accessing 
    %data: SD(2).depth(1,3)
    field1 = 'raw';         value1 = {SD1raw, SD2raw, SD3raw, ExtraSDraw};
    field2 = 'depth';       value2 = {SD1, SD2, SD3, ExtraSD(:,2:end)};
    field3 = 'glacier';     value3 = {SD1_Glacier, SD2_Glacier, SD3_Glacier, ExtraSD_Glacier};
    field4 = 'pattern';     value4 = {SD1_Pattern, SD2_Pattern, SD3_Pattern, ExtraSD_Pattern};
    field5 = 'person';      value5 = {SD1_Person, SD2_Person, SD3_Person, ExtraSD_Person};
    field6 = 'Q';           value6 = {SD1_Q, SD2_Q, SD3_Q, ExtraSD_Q};
    field7 = 'book';        value7 = {SD1_Book, SD2_Book, SD3_Book, ExtraSD_Book};
    field8 = 'comments';    value8 = {SD1comments, SD2comments, SD3comments, ExtraSDcomments};

    SD = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
        field5,value5,field6,value6,field7,value7,field8,value8);
            clear value* field* SD1* SD2* SD3* ExtraSD* gps_elev