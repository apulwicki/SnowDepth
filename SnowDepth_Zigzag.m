%% Working with zigzag data

%% Import data
run SnowDepth_Import.m %Imports snow depth and measurement location data
run SnowDensity_Import.m %Imports snow density values

%Compute distance of each measurement from its vertex
glacier_categories = categories(ZZ.glacier);
zone_categories = categories(ZZ.zone);
for k = 1:length(glacier_categories)
    glacier = glacier_categories(k);
    for j = 1:length(zone_categories)
        zone = zone_categories(j);
        indpre = intersect(find(ZZ.glacier==glacier),find(ZZ.zone==zone));
        for i = 1:8
            vertex = strcat('ZZ0',num2str(i));
            ind = intersect(indpre,find(ZZ.vertexlabel==vertex));
            if isempty(ind)
                continue
            else
                ZZ.depth(ind(1):ind(end),1) = cumsum(ZZ.depth(ind,2));
            end
        end
    end
end


ZZ_cord = [strcat(ZZ.text(:,1),'_',ZZ.text(:,2),'_',ZZ.text(:,3)),num2cell(ZZ.depth(:,1))];
    %use cell2mat to convert back to number

for j = 1:length(ZZ_cord)
    temp = strcmp(ZZ_cord(j,1),ZZ.vertexcoord(:,5));
    tempind = find(temp,1);
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
clear distinterp index i j k line temp tempind easting northing ind indpre ...
    zone glacier vertex zone_categories glacier_categories

%Adding data to ZZ structure
ZZ_cord = [ZZ_cord, num2cell(ZZ.depth(:,3:4))];
ZZ.depth = ZZ_cord; clear ZZ_cord 

%Adding measured density from SWE values
GZZindex = [1,150,318,482,674,835,987,1144,1289,1457,1620];
for i = 1:size(GZZindex,2)-1
    ZZ.depth(GZZindex(i):GZZindex(i+1)-1,7) = num2cell(cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,5))*cell2mat(zigzagSWE(i,2)));
end

%Selecting only good quality (Q=1) data
GZZutm = ZZ.depth(GZZindex,3:4);
ZZ.depth(cell2mat(ZZ.depth(:,6))==0,:)=[];
for i = 1:size(GZZindex,2)
    GZZindex(1,i) = find(cell2mat(GZZutm(i,1)) == cell2mat(ZZ.depth(:,3:4)),1);
end
clear GZZutm i k 
%% Plotting zigzag data

GZZlabel = ['G04 Z3A'; 'G04 Z2A'; 'G04 Z5B'; 'G02 Z5C'; 'G02 Z7A';'G02 Z3B'; 'G13 Z7C';'G13 Z4C'; 'G13 Z3B'; 'G13 Z5A'];
min = nanmin(cell2mat(ZZ.depth(:,7)));
max = nanmax(cell2mat(ZZ.depth(:,7)));

for i = 1:size(GZZlabel,1)
    x = cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,3));
    y = cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,4));
    z = cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,7));
    x2 = nanmax(x)-x;
    y2 = nanmax(y)-y;
    pointsize = 30;
    
    meandepth = mean(z);
    stddepth = std(z);
    
    figure(1)
    %subplot(2,2,i)
    scatter(x2, y2, pointsize, z,'filled');
    axis equal
    str = {strcat('mean= ', num2str(round(meandepth,1)),'cm SWE'), ...
        strcat('std= ', num2str(round(stddepth,1)),'cm SWE')};
    title(GZZlabel(i,:))
    dim = [.13 .5 .3 .3];
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
    xlabel('Distance (m)')
    ylabel('Distance (m)')
        caxis([min max])
        c = colorbar;
        %colormap(flipud(cool))
        c.Label.String = 'SWE (cm)';
   
   filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/',GZZlabel(i,:),'same_scale');
   %filename = strcat('/home/glaciology1/Documents/Data/Plots/same_scale',GZZlabel(i,:));
   print(filename,'-dpng')
   clf
%     if i == 3
%         min = nanmin(cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,5)));
%         max = nanmax(cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,5)));
%         caxis([min max])
%         c = colorbar;
%         colormap(flipud(cool))
%         c.Label.String = 'Snow depth (cm)';
%     end
    
end
clear max min x y x2 y2 z i k c pointsize stddepth meandepth str dim filename GZZlabel

%dlmwrite('/home/glaciology1/Documents/QGIS/Data/TestPoints.csv',closest, 'delimiter', ',', 'precision', 9); %Write matrix with new waypoints to csv file for QGIS

%% Variogram - zigzag

GZZlabel = ['G04 Z3A'; 'G04 Z2A'; 'G04 Z5B'; 'G02 Z5C'; 'G02 Z7A';'G02 Z3B'; 'G13 Z7C';'G13 Z4C'; 'G13 Z3B'; 'G13 Z5A'];
min = nanmin(cell2mat(ZZ.depth(:,7)));
max = nanmax(cell2mat(ZZ.depth(:,7)));

for i = 1:size(GZZlabel,1)
    x = cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,3));
    y = cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,4));
    z = cell2mat(ZZ.depth(GZZindex(i):GZZindex(i+1)-1,7));
    x2 = nanmax(x)-x;
    y2 = nanmax(y)-y;
    pointsize = 30;
    d = variogram([x2 y2],z,'plotit',false,'nrbins',100);
    

    %variogram fit
    figure(1)
%     subplot(2,1,1)
%         scatter(x2,y2,pointsize,z,'filled'); box on;
%         title(GZZlabel(i,:))
%     subplot(2,1,2)
        h=d.distance;
        gammaexp = d.val;
        a0 = 15; % initial value: range 
        c0 = 0.1; % initial value: sill 
        [a,c,n] = variogramfit(h,gammaexp,a0,c0,[],...
                               'solver','fminsearchbnd',...
                               'nugget',0,...
                               'plotit',true);
        title(GZZlabel(i,:))
        
   %filename = strcat('/home/glaciology1/Documents/Data/Plots/variogram',GZZlabel(i,:));
   filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/',GZZlabel(i,:),'variogram');
   print(filename,'-dpng')
   clf

end
clear x y x2 y2 z i k c pointsize filename GZZlabel

% figure(1)
%     subplot(2,2,1)
%     scatter(x2,y2,4,z,'filled'); box on;
%     ylabel('y'); xlabel('x')
%     %title(GZZ_lab(i,:))
%     subplot(2,2,2)
%     hist(z,20)
%     ylabel('frequency'); xlabel('z')
%     title('histogram of z-values')
%     subplot(2,2,3)
%     d = variogram([x2 y2],z,'plotit',true,'nrbins',100);
%     title('Isotropic variogram')
%     subplot(2,2,4)
%     d2 = variogram([x2 y2],z,'plotit',true,'nrbins',100,'anisotropy',true);
%     title('Anisotropic variogram')                   
%% Length of measured vs vertex spacing

section_sum = [1];
for i = 1:size(ZZ.depth,1)-1
    if ~cellfun(@isequal,ZZ.depth(i,1),ZZ.depth(i+1,1))
        section_sum = [section_sum ; i+1];
    end
end
for i = 1:size(section_sum,1)-1
    section_sum(i,2) = cell2mat(ZZ.depth(section_sum(i+1,1)-1,2));
end
section_sum(end,:) = [];
display(strcat('Section sum = ', num2str(mean(section_sum(:,2)))))

vertex_spacing = [];
for i = 2:8
    vertex_spacing(i,1) = EuclideanDistance(cell2mat(ZZ.vertexcoord(i,1)),cell2mat(ZZ.vertexcoord(i,2)),...
        cell2mat(ZZ.vertexcoord(i+1,1)),cell2mat(ZZ.vertexcoord(i+1,2)));
end
vertex_spacing(1,:) = [];
vertex_spacing(4,:) = [];
display(strcat('Vertex spacing = ',num2str(mean(vertex_spacing(:,1)))))

plot(section_sum(:,2),'o')
hold on
plot([0 size(section_sum,1)],[mean(vertex_spacing(:,1)) mean(vertex_spacing(:,1))])
    ylabel('Distance (m)')
    legend('Measured spacing','Vertex spacing', 'Location','best')

