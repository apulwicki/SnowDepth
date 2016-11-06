
%run Import_Topo.m

%% Maps of topographic params for each glacier
header  = fieldnames(topo_full.G4);
units   = {'', '(degrees)','(m a.s.l)','','','(degrees)',''};
glacier = {'G4','G2','G13'}; 
for r = 1:length(header)
figure(1)
    param = char(header(r));
    x_min = nanmin([topo_full.G4.(param)(:);topo_full.G2.(param)(:);topo_full.G13.(param)(:)]);
    x_max = nanmax([topo_full.G4.(param)(:);topo_full.G2.(param)(:);topo_full.G13.(param)(:)]);
    for i = 1:3
        name    = char(glacier(i)); 
        data    = topo_full.(name).(param);
        data(isnan(data)) = x_max;
            s(i) = subplot(1,3,i);
                imagesc(data); caxis([x_min max([x_max,max(data)])]);
                colordata = colormap; colordata(end,:) = [1 1 1]; colormap(colordata); % for making NaN values white
                axis square
            if i ==2; title(param); end            
            if i ==3; colorbar('location','eastoutside'); ylabel(c,units(r)); end
    end
            s1Pos = get(s(1),'position');
            s3Pos = get(s(3),'position');
            s3Pos(3:4) = s1Pos(3:4);
            set(s(3),'position',s3Pos);

        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = ['Map_',header{r}];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')

clf
end 

    clear i r name header glacier fig param s* x* data colordata
    
%% SWE at sampling locations

rig.G4 = shaperead('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G4.shp','UseGeoCoords', true);
rig.G2 = shaperead('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G2.shp','UseGeoCoords', true);
rig.G13 = shaperead('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G12.shp','UseGeoCoords', true);

pointsize = 13;

%%
for i = 1:3
    name    = char(glacier(i));
    s(i) = subplot(1,3,i);
        axesm utm; setm(gca,'zone','7v');  
        %h = getm(gca); setm(gca,'grid','on','meridianlabel','on','parallellabel','on')
        geoshow(rig.(name), 'FaceColor', [1 1 1]); hold on
        scatter(SWE(i).utm(:,1), SWE(i).utm(:,2), pointsize, SWE(i).swe,'filled'); 
            axis off; axis tight
            s(i).Box = 'off'; 
    if i ==3; c = colorbar('location','eastoutside'); ylabel(c,'SWE (m)'); end        
end
    linkprop(s,'Clim');
    s1Pos = get(s(1),'position');	s2Pos = get(s(2),'position');  s3Pos = get(s(3),'position');    
    s3Pos(3:4) = s2Pos(3:4);        set(s(3),'position',s3Pos);
    s1Pos(3:4) = s2Pos(3:4);        set(s(1),'position',s1Pos);

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = 'SWEmap';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
    
