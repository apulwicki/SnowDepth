
%run Import_Topo.m

%% Maps of topographic params for each glacier
header  = fieldnames(topo_full.G4);
pos_axis = [0 0 .3 1; 0.2 0 .3 1; 0.5 0 .3 1];

for r = 1:length(header)
figure(1)
    param = char(header(r));
    x_min = nanmin([topo_full.G4.(param)(:);topo_full.G2.(param)(:);topo_full.G13.(param)(:)]);
    x_max = nanmax([topo_full.G4.(param)(:);topo_full.G2.(param)(:);topo_full.G13.(param)(:)]);
    for i = 1:3
        name    = char(options.glacier(i)); 
        data    = topo_full.(name).(param);
        data(isnan(data)) = x_max;
            s(i) = subplot(1,3,i);
                imagesc(data); caxis([x_min max([x_max,max(data)])]); hold on
                colordata = colormap; colordata(end,:) = [1 1 1]; colormap(colordata); % for making NaN values white
                %!!!!!!plot(SWE(i).utm(:,1),SWE(i).utm(:,2),'k.');
                axis square
                s(i).Position   = pos_axis(i,:);
                s(i).Box        = 'off';  axis off
            if i ==2; title(param); end            
            if i ==3; c = colorbar('location','eastoutside');  ylabel(c,char(options.topoVars(r))); end 
    end
            s1Pos = get(s(1),'position');   s3Pos = get(s(3),'position');   
            s3Pos(3:4) = s1Pos(3:4);        set(s(3),'position',s3Pos);

        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = ['Map_',header{r}];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')

clf
end 

    clear i r name header glacier fig param s* x* data colordata
    
%% SWE at sampling locations

% rig.G4 = shaperead('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G4.shp','UseGeoCoords', true);
% rig.G2 = shaperead('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G2.shp','UseGeoCoords', true);
% rig.G13 = shaperead('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G12.shp','UseGeoCoords', true);

% rig.G4 = shaperead('/home/glaciology1/Documents/Data/GlacierShapeFiles/RIG_G4.shp','UseGeoCoords', true);
% rig.G2 = shaperead('/home/glaciology1/Documents/Data/GlacierShapeFiles/RIG_G2.shp','UseGeoCoords', true);
% rig.G13 = shaperead('/home/glaciology1/Documents/Data/GlacierShapeFiles/RIG_G12.shp','UseGeoCoords', true);

rig.G4 = csvread('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G4.csv');
rig.G2 = csvread('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G2.csv');
rig.G13 = csvread('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G13.csv');

pointsize = 13;
glacier = {'G4','G2','G13'};  

pos_axis = [0 0 .3 1; 0.2 0 .3 1; 0.5 0 .3 1];


Xlim = [0 max([rig.G4(:,1)-min(rig.G4(:,1));rig.G2(:,1)-min(rig.G2(:,1));rig.G13(:,1)-min(rig.G13(:,1))])];
Ylim = [0 max([rig.G4(:,2)-min(rig.G4(:,2));rig.G2(:,2)-min(rig.G2(:,2));rig.G13(:,2)-min(rig.G13(:,2))])];

for i = 1:2
    name    = char(glacier(i));
    s(i) = subplot(1,3,i);
            minE = min(rig.(name)(:,1));
            minN = min(rig.(name)(:,2));
        Eg = rig.(name)(:,1) - minE;
        Ng = rig.(name)(:,2) - minN;
        plot(Eg,Ng,'k'); hold on
        Es = SWE(i).utm(:,1) - minE;
        Ns = SWE(i).utm(:,2) - minN;
        scatter(Es,Ns , pointsize, SWE(i).swe,'filled');
            axis off;       axis equal;
            s(i).Box        = 'off'; 
            s(i).Position   = pos_axis(i,:);
        Xlim(2) = max([Xlim(2), max(Eg)]); Ylim(2) = max([Ylim(2), max(Ng)]); 
        axis([Xlim, Ylim]); caxis([0 1.2]);
end 

i = 3;
name    = char(glacier(i));
s(i) = subplot(1,3,i);
    minE = min(rig.(name)(:,1));
    minN = min(rig.(name)(:,2));
Eg = rig.(name)(:,1) - minE;
Ng = rig.(name)(:,2) - minN;
plot(Eg(1:304),Ng(1:304),'k'); hold on
plot(Eg(305:330),Ng(305:330),'k');
plot(Eg(331:356),Ng(331:356),'k');
plot(Eg(357:end),Ng(357:end),'k');
Es = SWE(i).utm(:,1) - minE;
Ns = SWE(i).utm(:,2) - minN;
scatter(Es,Ns , pointsize, SWE(i).swe,'filled');
    axis off;       axis equal;
    s(i).Box        = 'off'; 
    s(i).Position   = pos_axis(i,:);
Xlim(2) = max([Xlim(2), max(Eg)]); Ylim(2) = max([Ylim(2), max(Ng)]); 
axis([Xlim, Ylim]); caxis([0 1.2]);



if i ==3; c = colorbar('location','eastoutside'); 
    ylabel(c,'SWE (m)'); end   
        
        
        

    s2Pos = get(s(2),'position');   s3Pos = get(s(3),'position');    
    s3Pos(3:4) = s2Pos(3:4);        set(s(3),'position',s3Pos);

    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = 'SWEmap';
% print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
    

% scaleruler('units','km')
% setm(handlem('scaleruler'), ...
%     'MajorTick', 0:0.5:2,...
%     'MinorTick', 0,...
%     'TickDir', 'down', ...
%     'MajorTickLength', 0.25,...
%     'MinorTickLength', 0)


