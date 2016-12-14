
%run Import_Topo.m

%% Maps of topographic params for each glacier

header      = fieldnames(topo_full_ns.G4);

for r = 1:length(header)
    param = char(header(r));
    topoParam.G4 = topo_full_ns.G4.(param);
    topoParam.G2 = topo_full_ns.G2.(param);
    topoParam.G13 = topo_full_ns.G13.(param);
    
    PlotTopoParameter(topoParam,param, options.topoVarsUnits(r), SWE)
end 
    
%% SWE at sampling locations

% rig.G4 = csvread('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G4.csv');
% rig.G2 = csvread('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G2.csv');
% rig.G13 = csvread('/Users/Alexandra/Documents/SFU/Data/glaciershapefiles/RIG_G13.csv');

rig.G4 = csvread('/home/glaciology1/Documents/Data/GlacierShapeFiles/RIG_G4.csv');
rig.G2 = csvread('/home/glaciology1/Documents/Data/GlacierShapeFiles/RIG_G2.csv');
rig.G13 = csvread('/home/glaciology1/Documents/Data/GlacierShapeFiles/RIG_G13.csv');
rig.ELA = csvread('/home/glaciology1/Documents/Data/GlacierShapeFiles/ELA.csv');


pointsize = 13;

clf
pos_axis        = [0 0 .35 1; 0.19 0 .35 1; 0.51 0 .35 1];
flowarrow_x     = [0.07 0.08; 0.46 0.42; 0.74 0.73];
flowarrow_y     = [0.45 0.36; 0.23  0.31; 0.20 0.31];

for opt = 2:9

    run OPTIONS.m
    options.DensitySWE  = opt;
    options.ZZ          = 1; %include zigzags
    run MAIN

ELA_d = [1 7; 16 23; 8 15];
    
Xlim = [0 max([rig.G4(:,1)-min(rig.G4(:,1));rig.G2(:,1)-min(rig.G2(:,1));rig.G13(:,1)-min(rig.G13(:,1))])];
Ylim = [0 max([rig.G4(:,2)-min(rig.G4(:,2));rig.G2(:,2)-min(rig.G2(:,2));rig.G13(:,2)-min(rig.G13(:,2))])];
Zlim = [0 1.2]; %round(max([SWE(1).swe; SWE(2).swe; SWE(3).swe]),2)

for i = 1:2
    name    = char(options.glacier(i));
    s(i) = subplot(1,3,i);
       %Glacier outline
            minE = min(rig.(name)(:,1));
            minN = min(rig.(name)(:,2));
        Eg = rig.(name)(:,1) - minE;
        Ng = rig.(name)(:,2) - minN;
        plot(Eg,Ng,'k'); hold on
       %ELA
        Eela = rig.ELA(ELA_d(i,1):ELA_d(i,2),1) - minE;
        Nela = rig.ELA(ELA_d(i,1):ELA_d(i,2),2) - minN;
        plot(Eela,Nela ,'k--')
       %SWE
        Es = SWE(i).utm(:,1) - minE;
        Ns = SWE(i).utm(:,2) - minN;
        scatter(Es,Ns , pointsize, SWE(i).swe,'filled');
            axis off;       axis equal;
            s(i).Box        = 'off'; 
            s(i).Position   = pos_axis(i,:);
        Xlim(2) = max([Xlim(2), max(Eg)]); Ylim(2) = max([Ylim(2), max(Ng)]); 
        axis([Xlim, Ylim]); caxis([0 1.2]);
        annotation('textarrow',flowarrow_x(i,:),flowarrow_y(i,:),'String','')
%Scale bar
if i ==1; scalebar('ScaleLength', 2000, 'Unit', 'm', 'Location','northwest'); end

end 

i = 3;
    name    = char(options.glacier(i));
    s(i) = subplot(1,3,i);
    %Glacier outline
        minE = min(rig.(name)(:,1));
        minN = min(rig.(name)(:,2));
    Eg = rig.(name)(:,1) - minE;
    Ng = rig.(name)(:,2) - minN;
    plot(Eg(1:304),Ng(1:304),'k'); hold on
    plot(Eg(305:330),Ng(305:330),'k');
    plot(Eg(331:356),Ng(331:356),'k');
    plot(Eg(357:end),Ng(357:end),'k');
       %ELA
        Eela = rig.ELA(ELA_d(i,1):ELA_d(i,2),1) - minE;
        Nela = rig.ELA(ELA_d(i,1):ELA_d(i,2),2) - minN;
        plot(Eela,Nela ,'k--')
    %SWE
    Es = SWE(i).utm(:,1) - minE;
    Ns = SWE(i).utm(:,2) - minN;
    s3 = scatter(Es,Ns , pointsize, SWE(i).swe,'filled');
        axis off;       axis equal;
        s(i).Box        = 'off'; 
        s(i).Position   = pos_axis(i,:);
        annotation('textarrow',flowarrow_x(i,:),flowarrow_y(i,:),'String','')
    Xlim(2) = max([Xlim(2), max(Eg)]); Ylim(2) = max([Ylim(2), max(Ng)]); 
    axis([Xlim, Ylim]); caxis(Zlim);

    if i ==3; c = colorbar('location','eastoutside');
        ylabel(c,'SWE (m)'); end   
  

     % North arrow
    Narrow = imread('Narrow.jpg');
    a = axes('position',[-0.02,0.65,0.12,0.12]); 
    Nshow = imshow(Narrow, parula);   
    colormap(a,gray)

       

    s2Pos = get(s(2),'position');   s3Pos = get(s(3),'position');    
    s3Pos(3:4) = s2Pos(3:4);        set(s(3),'position',s3Pos);

    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];
filename = ['SWEmap_opt',num2str(opt)];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
    
end

%% Modelled and observed SWE
modelled = sweMLR;
type     = 'MLR';

pos_axis    = [0 0 .4 1; 0.19 0 .4 1; 0.50 0 .4 1];
opt = 8;

    x_min = nanmin([modelled(opt).G4(:);modelled(opt).G2(:);modelled(opt).G13(:)]);
    x_max = nanmax([modelled(opt).G4(:);modelled(opt).G2(:);modelled(opt).G13(:)]);

figure(1)      
for i = 1:3
    name    = char(options.glacier(i)); 
    G13size = size(modelled(opt).G13);
    Gsize   = size(modelled(opt).(name));
    data    = [ nan(G13size(1,1)-Gsize(1,1),Gsize(1,2)); modelled(opt).(name)];
    data    = [data, nan(G13size(1,1), G13size(1,2)-Gsize(1,2))];

    s(i) = subplot(1,3,i);
        h = imagesc(data); hold on
            caxis([x_min max([x_max,max(data)])]); 
            set(h,'alphadata',~isnan(data))
            E = (SWE(i).utm(:,1)-min(SWE(i).utm(:,1)))/40;
            Na = SWE(i).utm(:,2)-min(SWE(i).utm(:,2)); N = (max(Na)-Na)/40;                
            if      i==1; E = E+25; N = N+78;
            elseif  i==2; E = E+10; N = N+47;
            elseif  i==3; E = E+17; N = N+21; 
            end
        %SWE
        Ns = (SWE(i).utm(:,2) - min(rig.(name)(:,2)))/40-14;  Ns = 140-Ns;
        Es = (SWE(i).utm(:,1) - min(rig.(name)(:,1)))/40; 
        plot(Es,Ns,'k.')
        scatter(Es,Ns , 13, SWE(i).swe,'filled');            
            axis equal
            s(i).Position   = pos_axis(i,:);
            s(i).Box        = 'off';  axis off
            %Scale bar
            if      i ==1; scalebar('ScaleLength', 2000, 'Unit', 'm', 'Location','northwest');  
            elseif  i ==3; c = colorbar('location','eastoutside');  ylabel(c,'Modelled SWE (m w.e.)'); 
            end 
end
            s1Pos = get(s(1),'position');   s3Pos = get(s(3),'position');   
            s3Pos(3:4) = s1Pos(3:4);        set(s(3),'position',s3Pos);

     % North arrow
    Narrow = imread('Narrow.jpg');
    a = axes('position',[-0.02,0.65,0.12,0.12]); 
    Nshow = imshow(Narrow, parula);   
    colormap(a,gray)    

    set(c,'Position',[0.84 0.2 0.036 0.6])

        
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];
filename = [type,'map_Modelled_Observed',char(opt)];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')

    clear i r name header glacier fig param s x* data colordata filename Es a c N*
    clear data E g G13size Gsize h modelled N Na name opt pos_axis t type s1Pos s3Pos
    
%% Modelled SWE Difference
modelled = sweBMS;

    %Differencing modelled SWE
for g = 1:3
    glacier = char(options.glacier(g));
    
    for i = 2:9
    stackSWE.(glacier)(:,:,i-1)   = modelled(i).(glacier);
    end

minSWE.(glacier)  = nanmin(stackSWE.(glacier),[],3);
maxSWE.(glacier)  = nanmax(stackSWE.(glacier),[],3);

diffSWE.(glacier) = maxSWE.(glacier)-minSWE.(glacier);
end
    
    %Maps
pos_axis    = [0 0 .4 1; 0.19 0 .4 1; 0.50 0 .4 1];
    x_min = nanmin([diffSWE.G4(:);diffSWE.G2(:);diffSWE.G13(:)]);
    x_max = nanmax([diffSWE.G4(:);diffSWE.G2(:);diffSWE.G13(:)]);
        
for i = 1:3
    name    = char(options.glacier(i)); 
    G13size = size(diffSWE.G13);
    Gsize   = size(diffSWE.(name));
    data    = [ nan(G13size(1,1)-Gsize(1,1),Gsize(1,2)); diffSWE.(name)];
    data    = [data, nan(G13size(1,1), G13size(1,2)-Gsize(1,2))];

    s(i) = subplot(1,3,i);
            h = imagesc(data); hold on
            caxis([x_min max([x_max,max(data)])]); 
            set(h,'alphadata',~isnan(data))
            axis equal
            s(i).Position   = pos_axis(i,:);
            s(i).Box        = 'off';  axis off
        %Scale bar
        if      i ==1; scalebar('ScaleLength', 2000, 'Unit', 'm', 'Location','northwest');  
        %elseif  i ==2; title(param);             
        elseif  i ==3; c = colorbar('location','eastoutside');  ylabel('Difference in modelled SWE values (m w.e)'); 
        end 
end
        s1Pos = get(s(1),'position');   s3Pos = get(s(3),'position');   
        s3Pos(3:4) = s1Pos(3:4);        set(s(3),'position',s3Pos);

     % North arrow
    Narrow = imread('Narrow.jpg');
    a = axes('position',[-0.02,0.65,0.12,0.12]); 
    Nshow = imshow(Narrow, parula);   
    colormap(a,gray)    

    set(c,'Position',[0.84 0.2 0.036 0.6])

        
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];
filename = 'SWEdifferenceMap';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
 


