%%%%%%%%%%%%%%%%%
% Importing Snow Depth data
%%%%%%%%%%%%%%%%%

%% Importing data and categorizing

%Import data from excel sheet
[SD1, SD1text, SD1raw] = xlsread('Field Data.xls','SD#1'); %Import values, text, and raw data from first sheet
SD1header = SD1text(1,:); %Headers for first sheet are in row 1
SD1text = SD1text(2:end,:); %Remaining text in tsheet 1
[SD2, SD2text, SD2raw] = xlsread('Field Data.xls','SD#2');
SD2header = SD2text(1,:);
SD2text = SD2text(2:end,:);
[SD3, SD3text, SD3raw] = xlsread('Field Data.xls','SD#3'); 
SD3header = SD3text(1,:);
SD3text = SD3text(2:end,:);
[ZZ, ZZtext, ZZraw] = xlsread('Field Data.xls','ZigZag'); 
ZZheader = ZZtext(1,:);
ZZtext = ZZtext(2:end,:);

%Categorize quality and situational data
SD1_Q1 = categorical(SD1(:,3));
SD1_Q2 = categorical(SD1(:,5));
SD1_Q3 = categorical(SD1(:,7));
SD1_Q4 = categorical(SD1(:,9));
SD1_Book = categorical(SD1text(:,11));
SD1_Glacier = categorical(SD1text(:,12));
SD1_Person = categorical(SD1text(:,13));
SD1_Pattern = categorical(SD1text(:,14));
SD1_Date = categorical(SD1text(:,15));

SD2_Q1 = categorical(SD2(:,3));
SD2_Q2 = categorical(SD2(:,5));
SD2_Q3 = categorical(SD2(:,7));
SD2_Q4 = categorical(SD2(:,9));
SD2_Book = categorical(SD2text(:,11));
SD2_Glacier = categorical(SD2text(:,12));
SD2_Person = categorical(SD2text(:,13));
SD2_Pattern = categorical(SD2text(:,14));
SD2_Date = categorical(SD2text(:,15));

SD3_Q1 = categorical(SD3(:,3));
SD3_Q2 = categorical(SD3(:,5));
SD3_Q3 = categorical(SD3(:,7));
SD3_Q4 = categorical(SD3(:,9));
SD3_Book = categorical(SD3text(:,11));
SD3_Glacier = categorical(SD3text(:,12));
SD3_Person = categorical(SD3text(:,13));
SD3_Pattern = categorical(SD3text(:,14));
SD3_Date = categorical(SD3text(:,15));

ZZ_Glacier = categorical(ZZtext(:,1));
ZZ_Zone = categorical(ZZtext(:,2));
ZZ_Vertex = categorical(ZZtext(:,3));
ZZ_Q = categorical(ZZtext(:,6));
ZZ_Person = categorical(ZZtext(:,8));
ZZ_Date = categorical(ZZtext(:,7));

%% 
