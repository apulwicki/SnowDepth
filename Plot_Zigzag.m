%% Plotting zigzag data

GZZlabel = ['G04 Z3A'; 'G04 Z2A'; 'G04 Z5B'; 'G02 Z5C'; 'G02 Z7A';'G02 Z3B'; 'G13 Z7C';'G13 Z4C'; 'G13 Z3B'; 'G13 Z5A'];
min1 = nanmin(cell2mat(ZZ.depth(:,5)));
max1 = nanmax(cell2mat(ZZ.depth(:,5)));

for i = 1:size(GZZlabel,1)
    x = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,3));
    y = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,4));
    z = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,5));
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
        caxis([min1 max1])
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
clear max1 min1 x y x2 y2 z i k c pointsize stddepth meandepth str dim filename GZZlabel

%dlmwrite('/home/glaciology1/Documents/QGIS/Data/TestPoints.csv',closest, 'delimiter', ',', 'precision', 9); %Write matrix with new waypoints to csv file for QGIS

%% Variogram - zigzag

GZZlabel = ['G04 Z3A'; 'G04 Z2A'; 'G04 Z5B'; 'G02 Z5C'; 'G02 Z7A';'G02 Z3B'; 'G13 Z7C';'G13 Z4C'; 'G13 Z3B'; 'G13 Z5A'];
min1 = nanmin(cell2mat(ZZ.depth(:,5)));
max1 = nanmax(cell2mat(ZZ.depth(:,5)));

for i = 1:size(GZZlabel,1)
    x = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,3));
    y = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,4));
    z = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,5));
    x2 = nanmax(x)-x;
    y2 = nanmax(y)-y;
    d = variogram([x2 y2],z,'plotit',false,'nrbins',100, 'maxdist', 40);
    

    %variogram fit
    figure(1)
%     subplot(2,1,1)
%         pointsize = 30;
%         scatter(x2,y2,pointsize,z,'filled'); box on;
%         title(GZZlabel(i,:))
     subplot(3,1,1:2)
        h=d.distance;
        gammaexp = d.val;
        a0 = 15; % initial value: range 
        c0 = 0.1; % initial value: sill 
        [a,c,n,S] = variogramfit(h,gammaexp,a0,c0,[],...
                               'solver','fminsearchbnd',...
                               'nugget',0,'plotit',true,...
                               'model','spherical');
        str = ['R^2 = ',num2str(round(S.Rs,3))];
        t = annotation('textbox',[.17 .7 .2 .2],'string',str,'FitBoxToText','on');
        s = t.FontSize;
        t.FontSize = 14;
        title(GZZlabel(i,:))
        
        
        subplot(3,1,3)
        lag = 5;
        group_dist = (0:lag:round(max(d.distance),-1))';
        group_num = zeros(size(group_dist,1)-1,1);
        for j = 1:size(group_dist,1)-1
            index = intersect(find(d.distance>group_dist(j,1)), find(d.distance<group_dist(j+1,1)));
            group_num(j,1) = sum(d.num(index,1));
        end
        bar(mean([group_dist(1:end-1), group_dist(2:end)],2),group_num,'BarWidth', 1)
        xlabel('Lag'); ylabel('# Pairs');
        
        axes('Position',[.71 .46 .15 .13])
        box on
        plot(d.distance,d.num)
        ylabel('# pairs'); xlabel('lag');
       
   %filename = strcat('/home/glaciology1/Documents/Data/Plots/variogram',GZZlabel(i,:));
    filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/',GZZlabel(i,:),'variogram');
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 9];
    print(filename,'-dpng','-r0')
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

