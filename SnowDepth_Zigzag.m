%% Working with zigzag data

%Compute distance of each measurement from its vertex
glacier_categories = categories(ZZ_Glacier);
zone_categories = categories(ZZ_Zone);
for k = 1:length(glacier_categories)
    glacier = glacier_categories(k);
    for j = 1:length(zone_categories)
        zone = zone_categories(j);
        indpre = intersect(find(ZZ_Glacier==glacier),find(ZZ_Zone==zone));
        for i = 1:8
            vertex = strcat('ZZ0',num2str(i));
            ind = intersect(indpre,find(ZZ_Vertex==vertex));
            if isempty(ind)
                continue
            else
                ZZ(ind(1):ind(end),1) = cumsum(ZZ(ind,2));
            end
        end
    end
end


ZZ_cord = [strcat(ZZtext(:,1),'_',ZZtext(:,2),'_',ZZtext(:,3)),num2cell(ZZ(:,1))];
    %use cell2mat to convert back to number

for j = 1:length(ZZ_cord)
    temp = strcmp(ZZ_cord(j,1),Vertex_cord(:,5));
    tempind = find(temp,1);
    easting = [cell2mat(Vertex_cord(tempind,1)), cell2mat(Vertex_cord(tempind+1,1))];
    northing = [cell2mat(Vertex_cord(tempind,2)), cell2mat(Vertex_cord(tempind+1,2))];

    if cell2mat(strfind(Vertex_cord(tempind,5),'8'))==12
        easting = [cell2mat(Vertex_cord(tempind,1)), cell2mat(Vertex_cord(tempind-7,1))];
        northing = [cell2mat(Vertex_cord(tempind,2)), cell2mat(Vertex_cord(tempind-7,2))];
    end
     % Create a series of points between two GPS locations 
    line = linspaceNDim([easting(1,1),northing(1,1)],[easting(1,2),northing(1,2)],1000);

    % Find the distance between created points and the original GPS location
    distinterp = zeros(length(line),1);
    for i = 1:length(line)
        distinterp(i,1) = EuclideanDistance(line(1,i),line(2,i),easting(1,1),northing(1,1));
    end

    [~, index] = nanmin(abs(distinterp-cell2mat(ZZ_cord(j,2))));
    ZZ_cord(j,3:4) = num2cell(line(:,index)');

end
clear distinterp index i j line temp tempind easting northing ind indpre zone glacier vertex

ZZ_cord = [ZZ_cord, num2cell(ZZ(:,3))];

%% Plotting zigzag data
GZZ = [1,153,321,489,685,848,1004,1165,1314,1485,1653];
GZZ_lab = ['G04 Z3A'; 'G04 Z2A'; 'G04 Z5B'; 'G02 Z5C'; 'G02 Z7A';'G02 Z3B'; 'G13 Z7C';'G13 Z4C'; 'G13 Z3B'; 'G13 Z5A'];
for i = 1:3
    x = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,3));
    y = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,4));
    z = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,5));
    x2 = nanmax(x)-x;
    y2 = nanmax(y)-y;
    pointsize = 20;
    
    figure(1)
    subplot(2,2,i)
    scatter(x2, y2, pointsize, z,'filled');
    axis equal
    title(GZZ_lab(i,:))
    xlabel('Distance (m)')
    ylabel('Distance (m)')
    if i == 3
        min = nanmin(cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,5)));
        max = nanmax(cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,5)));
        caxis([min max])
        c = colorbar;
        colormap(flipud(cool))
        c.Label.String = 'Snow depth (cm)';
    end
    
end

clear max min x y x2 y2 z i k c

i=10;
x = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,3));
y = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,4));
z = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,5));
x2 = nanmax(x)-x;
y2 = nanmax(y)-y;

figure(1)
    scatter(x2, y2, pointsize, z,'filled');
    axis equal
    title(GZZ_lab(i,:))
    xlabel('Distance (m)')
    ylabel('Distance (m)')
        c = colorbar;
        colormap(flipud(cool))
        c.Label.String = 'Snow depth (cm)';

%dlmwrite('/home/glaciology1/Documents/QGIS/Data/TestPoints.csv',closest, 'delimiter', ',', 'precision', 9); %Write matrix with new waypoints to csv file for QGIS

%% Variogram - zigzag

i=10;
x = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,3));
y = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,4));
z = cell2mat(ZZ_cord(GZZ(i):GZZ(i+1)-1,5));
x2 = nanmax(x)-x;
y2 = nanmax(y)-y;

figure(1)
subplot(2,2,1)
scatter(x2,y2,4,z,'filled'); box on;
ylabel('y'); xlabel('x')
title(GZZ_lab(i,:))
subplot(2,2,2)
hist(z,20)
ylabel('frequency'); xlabel('z')
title('histogram of z-values')
subplot(2,2,3)
d = variogram([x2 y2],z,'plotit',true,'nrbins',100);
title('Isotropic variogram')
subplot(2,2,4)
d2 = variogram([x2 y2],z,'plotit',true,'nrbins',100,'anisotropy',true);
title('Anisotropic variogram')

%variogram fit
figure(2)
h=d.distance;
gammaexp = d.val;
a0 = 15; % initial value: range 
c0 = 0.1; % initial value: sill 
[a,c,n] = variogramfit(h,gammaexp,a0,c0,[],...
                       'solver','fminsearchbnd',...
                       'nugget',0,...
                       'plotit',true);
                   