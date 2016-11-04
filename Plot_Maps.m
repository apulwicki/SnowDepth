
%run Import_Topo.m

%% Maps of topographic params for each glacier
header = fieldnames(topo_full.G4);
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
            if i ==3; 
                colorbar('location','eastoutside'); end
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

pointsize = 13;
    x_min = nanmin([SWE(1).swe;SWE(2).swe;SWE(3).swe]);
    x_max = nanmax([SWE(1).swe;SWE(2).swe;SWE(3).swe]);

for i = 1:3
    E = SWE(i).utm(:,1)-min(SWE(i).utm(:,1));
    N = SWE(i).utm(:,2)-min(SWE(i).utm(:,2));
    subplot(1,3,i)
        s(i) = scatter(E, N, pointsize, SWE(i).swe,'filled');
            caxis([x_min x_max]);
        	axis square
    if i ==3; 
            colorbar('location','eastoutside'); end        
end
    s1Pos = get(s(1),'position');
    s3Pos = get(s(3),'position');
    s3Pos(3:4) = s1Pos(3:4);
    set(s(3),'position',s3Pos);


%!!!
S = shaperead(filename);

    
    