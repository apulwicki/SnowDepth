% Importing Zigzag Snow Depth Data
%       This script imports snowdepth data from the field data file for
%       transects, zigzags, and SWE. It also categorizes the descriptive
%       part of the data. Then it creates a structured array for zigzag data
%
%       Inputs:         Field Data ('FieldDataRevisedAP.xlsx')
%                       Zigzag corners 
%       Outputs:        Zigzag structure (ZZ)

%       Alexandra Pulwicki  Created: July 2016
%                           Updated: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing data

%Import data from excel sheet (change path based on computer)
    [ZZ, ZZtext, ZZraw] = xlsread('FieldDataRevisedAP.xlsx','ZigZag','A1:H1621'); 
    ZZtext(1,:) = [];

%Import vertex waypoints
    [~,~,Vertex_cord] = xlsread('zigzag_corners_utm.xls');

%% Categorizing data

%Categorize quality and situational data
    %Converting the text to categories means that it can be used to group
    %data for plotting and stats analysis in the future     

%Zigzag data
    ZZ_Glacier = categorical(ZZtext(:,1));                      %Glacier (G13, G02, G04)
    ZZ_Zone = categorical(ZZtext(:,2));                         %Zone label
    ZZ_Vertex = categorical(ZZtext(:,3));                       %Reference vertex label
    ZZ_Q = categorical(ZZ(:,3));                                %Quality of data (1 for good, 0 for question mark)
    ZZ_Person = categorical(ZZtext(:,8));                       %Person that took the measurement (AP, CA, GF, AC)
    ZZ_Date = categorical(ZZtext(:,7));                         %Date meaurement was taken
    C = cell(size(ZZ_Glacier));                                 
    C(:) = {'ZZ'};
    ZZ_Book = categorical(C);                                   %Book (all ZZ)
ZZ = [zeros(length(ZZ),1),ZZ];                              %[zeros(will be WP#), distance data, depth, quality]

clear C
%% Creating the structure for zigzag snowdepth data -> ZZ
    %Field values are names for columns, values are the cell arrays within the
    %structure. Only one row. Example of accessing data: ZZ.depth(3,1)
    field1 = 'raw';         value1 = {ZZraw};
    field2 = 'depth';       value2 = {ZZ};
    field3 = 'vertexlabel'; value3 = {ZZ_Vertex};
    field4 = 'vertexcoord'; value4 = {Vertex_cord};
    field5 = 'zone';        value5 = {ZZ_Zone};
    field6 = 'glacier';     value6 = {ZZ_Glacier};
    field7 = 'person';      value7 = {ZZ_Person};
    field8 = 'Q';           value8 = {ZZ_Q};
    field9 = 'book';        value9 = {ZZ_Book};
    field10 = 'text';       value10 = {ZZtext};

    ZZ = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
        field5,value5,field6,value6,field7,value7,field8,value8,field9,value9,field10,value10);
            clear value* field* ZZ_* ZZraw ZZtext Vertex_cord
            
%% Compute distance of each measurement from its vertex
    glacier_categories = categories(ZZ.glacier); %get names of glaciers
    zone_categories = categories(ZZ.zone); %get names of zones

% Get distance of each point from it's reference vertex (cumulative sum of
% distance for each line)
    for k = 1:length(glacier_categories) %for all glaciers
        glacier = glacier_categories(k); %set glacier name
        for j = 1:length(zone_categories) %for all zones
            zone = zone_categories(j); %set zone name
            indpre = intersect(find(ZZ.glacier==glacier),find(ZZ.zone==zone)); %get all points with desired glacier and zone 
            for i = 1:8 %there were 8 vertex points for each zigzag
                vertex = strcat('ZZ0',num2str(i)); %name of the vertex
                ind = intersect(indpre,find(ZZ.vertexlabel==vertex)); %find all points that started at the vertex
                if isempty(ind)
                    continue
                else
                    ZZ.depth(ind(1):ind(end),1) = cumsum(ZZ.depth(ind,2)); %get the distance of each point from its refernec vertex
                end
            end
        end
    end

    % Reference vertex and the distance from it for each point
    ZZ_cord = [strcat(ZZ.text(:,1),'_',ZZ.text(:,2),'_',ZZ.text(:,3)),num2cell(ZZ.depth(:,1))];
        %use cell2mat to convert back to number
if optionsZ.location == 1 
    for j = 1:length(ZZ_cord)
        tempind = find(strcmp(ZZ_cord(j,1),ZZ.vertexcoord(:,5)),1);
        easting = [cell2mat(ZZ.vertexcoord(tempind,1)), cell2mat(ZZ.vertexcoord(tempind+1,1))];
        northing = [cell2mat(ZZ.vertexcoord(tempind,2)), cell2mat(ZZ.vertexcoord(tempind+1,2))];

        if cell2mat(strfind(ZZ.vertexcoord(tempind,5),'8'))==12
            easting = [cell2mat(ZZ.vertexcoord(tempind,1)), cell2mat(ZZ.vertexcoord(tempind-7,1))];
            northing = [cell2mat(ZZ.vertexcoord(tempind,2)), cell2mat(ZZ.vertexcoord(tempind-7,2))];
        end

        x = easting(1,1) + cell2mat(ZZ_cord(j,2))*(easting(1,2)-easting(1,1))/...
            EuclideanDistance(easting(1,1),northing(1,1),easting(1,2),northing(1,2));
        y = northing(1,1) + cell2mat(ZZ_cord(j,2))*(northing(1,2)-northing(1,1))/...
            EuclideanDistance(easting(1,1),northing(1,1),easting(1,2),northing(1,2));

        ZZ_cord(j,3:4) = [num2cell(x), num2cell(y)];

    end
    
elseif optionsZ.location == 2
    use_vertex = ['G04_Z3A_ZZ01';'G04_Z3A_ZZ05';'G04_Z2A_ZZ01';'G04_Z2A_ZZ05';...
                   'G04_Z5B_ZZ01'; 'G04_Z5B_ZZ05';'G02_Z5C_ZZ08';'G02_Z7A_ZZ01';...
                   'G02_Z7A_ZZ08';'G02_Z7A_ZZ04';'G02_Z3B_ZZ03';'G02_Z3B_ZZ07';...
                   'G13_Z7C_ZZ02';'G13_Z7C_ZZ06';'G13_Z4C_ZZ08';'G13_Z4C_ZZ04';...
                   'G13_Z3B_ZZ04';'G13_Z3B_ZZ08';'G13_Z5A_ZZ08';'G13_Z5A_ZZ04'];
    
    for i = 1:length(ZZ_cord)
       if ismember(ZZ_cord(i,1),use_vertex)
            tempind = find(strcmp(ZZ_cord(i,1),ZZ.vertexcoord(:,5)),1);
            easting = [cell2mat(ZZ.vertexcoord(tempind,1)), cell2mat(ZZ.vertexcoord(tempind+1,1))];
            northing = [cell2mat(ZZ.vertexcoord(tempind,2)), cell2mat(ZZ.vertexcoord(tempind+1,2))];
            if cell2mat(strfind(ZZ.vertexcoord(tempind,5),'8'))==12
                easting = [cell2mat(ZZ.vertexcoord(tempind,1)), cell2mat(ZZ.vertexcoord(tempind-7,1))];
                northing = [cell2mat(ZZ.vertexcoord(tempind,2)), cell2mat(ZZ.vertexcoord(tempind-7,2))];
            end
       elseif ~ismember(ZZ_cord(i,1),ZZ_cord(i-1,1))            
           for j = 1:length(ZZ_cord)
               if isempty(cell2mat(ZZ_cord(j,3)))
                   tempind = find(strcmp(ZZ_cord(i,1),ZZ.vertexcoord(:,5)),1);
                   easting = [cell2mat(ZZ_cord(j-1,3)), cell2mat(ZZ.vertexcoord(tempind+1,1))];
                   northing = [cell2mat(ZZ_cord(j-1,4)), cell2mat(ZZ.vertexcoord(tempind+1,2))];
                   if cell2mat(strfind(ZZ.vertexcoord(tempind,5),'8'))==12
                        easting = [cell2mat(ZZ_cord(j-1,3)), cell2mat(ZZ.vertexcoord(tempind-7,1))];
                        northing = [cell2mat(ZZ_cord(j-1,4)), cell2mat(ZZ.vertexcoord(tempind-7,2))];
                   end
                   break
               end
           end
       end
       
        x = easting(1,1) + cell2mat(ZZ_cord(i,2))*(easting(1,2)-easting(1,1))/...
                EuclideanDistance(easting(1,1),northing(1,1),easting(1,2),northing(1,2));
        y = northing(1,1) + cell2mat(ZZ_cord(i,2))*(northing(1,2)-northing(1,1))/...
                EuclideanDistance(easting(1,1),northing(1,1),easting(1,2),northing(1,2));

        ZZ_cord(i,3:4) = [num2cell(x), num2cell(y)];   
       
    end
end    
%     scatter(cell2mat(ZZ_cord(:,3)),cell2mat(ZZ_cord(:,4))) %834
%     hold on
%     scatter(cell2mat(ZZ_cord(:,3)),cell2mat(ZZ_cord(:,4)), 'filled')
% %     hold on 
%      scatter(easting(1,:),northing(1,:),'filled')
%      hold on
%          scatter(cell2mat(ZZ_cord(699:726,3)),cell2mat(ZZ_cord(699:726,4)),'filled') %834

%%     
%     for j = 1:length(ZZ_cord)
%         tempind = find(strcmp(ZZ_cord(j,1),ZZ.vertexcoord(:,5)),1);
%         easting = [cell2mat(ZZ.vertexcoord(tempind,1)), cell2mat(ZZ.vertexcoord(tempind+1,1))];
%         northing = [cell2mat(ZZ.vertexcoord(tempind,2)), cell2mat(ZZ.vertexcoord(tempind+1,2))];
% 
%         if cell2mat(strfind(ZZ.vertexcoord(tempind,5),'8'))==12
%             easting = [cell2mat(ZZ.vertexcoord(tempind,1)), cell2mat(ZZ.vertexcoord(tempind-7,1))];
%             northing = [cell2mat(ZZ.vertexcoord(tempind,2)), cell2mat(ZZ.vertexcoord(tempind-7,2))];
%         end
% 
%         x = easting(1,1) + cell2mat(ZZ_cord(j,2))*(easting(1,2)-easting(1,1))/...
%             EuclideanDistance(easting(1,1),northing(1,1),easting(1,2),northing(1,2));
%         y = northing(1,1) + cell2mat(ZZ_cord(j,2))*(northing(1,2)-northing(1,1))/...
%             EuclideanDistance(easting(1,1),northing(1,1),easting(1,2),northing(1,2));
% 
%         ZZ_cord(j,3:4) = [num2cell(x), num2cell(y)];
% 
%     end
    
%end
clear distinterp index i j k line temp tempind easting northing ind indpre ...
    zone glacier vertex zone_categories glacier_categories

%Adding data to ZZ structure
ZZ_cord = [ZZ_cord, num2cell(ZZ.depth(:,3:4))];
ZZ.depth = ZZ_cord; clear ZZ_cord 

%Adding measured density from SWE values
if optionsZ.z == 2
    GZZindex = [1,150,318,482,674,835,987,1144,1289,1457,1620];
    for i = 1:size(GZZindex,2)-1
        ZZ.depth(GZZindex(i):GZZindex(i+1)-1,5) = num2cell(cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,5))*cell2mat(Density.zigzagtube(i,2)));
    end
end
    
%Selecting only good quality (Q=1) data
GZZutm = ZZ.depth(GZZindex,3:4);
ZZ.depth(cell2mat(ZZ.depth(:,6))==0,:)=[];
for i = 1:size(GZZindex,2)
    GZZindex(1,i) = find(cell2mat(GZZutm(i,1)) == cell2mat(ZZ.depth(:,3:4)),1);
end
ZZ.depth(:,6) = [];
clear GZZutm i k x y