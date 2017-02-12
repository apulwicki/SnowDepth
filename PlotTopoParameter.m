function [ ] = PlotTopoParameter( topoParam, paramName, cLabel, SWE, sweDOTS)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global options
    rig = topoParam.rig;
    ELA_d = [1 7; 16 23; 8 15];
clf
colormap('default')

%get colour min max
     %Manual colour range = 1
   fixed_color = 0;
if fixed_color == 0 
    x_min   = nanmin([topoParam.G4(:);topoParam.G2(:);topoParam.G13(:)]);
    x_max   = nanmax([topoParam.G4(:);topoParam.G2(:);topoParam.G13(:)]);
    minSWE  = nanmin([SWE(1).swe(:);SWE(2).swe(:);SWE(3).swe(:)]);
    maxSWE  = nanmax([SWE(1).swe(:);SWE(2).swe(:);SWE(3).swe(:)]);
elseif fixed_color == 1
    x_min = 0;      minSWE  = x_min;
    x_max = 1.4;    maxSWE  = x_max;
end

G13size = size(topoParam.G13);

width  = 0.41;
height = 2.3;
ystart = -0.64;
pos_axis    = [0        ystart   width     height; 
               0.15     ystart   width     height; 
               0.48     ystart   width     height];

%phasemap needs
if  strcmp(paramName, 'banana') %None need circular phase map, the option remains though
    Cmap = phasemap;
else
    Cmap = parula;
end
 
for i = 1:3
        name    = char(options.glacier(i)); 
        Gsize   = size(topoParam.(name));
        data    = [ nan(G13size(1,1)-Gsize(1,1),Gsize(1,2)); topoParam.(name)];
        data    = [data, nan(G13size(1,1), G13size(1,2)-Gsize(1,2))];
        
        s = subplot(1,3,i);
                if ~all(isnan(topoParam.(name)(:)))
                h = imagesc(data); hold on
                set(h,'alphadata',~isnan(data)); end
                
                if all(isnan(topoParam.(name)(:))) || strcmp(sweDOTS,'symmetric')
                        minE = min(rig.(name)(:,1));
                        minN = min(rig.(name)(:,2));
                    Eg = (rig.(name)(:,1) - minE)/40;
                    Ng = (rig.(name)(:,2) - minN)/40;  Ng = (max(Ng)-Ng);
                    Ng = Ng + (G13size(1,1)-max(Ng));
                    if i ==3 %Deals with G13 multiple parts
                        plot(Eg(1:304),Ng(1:304),'k'); hold on
                        plot(Eg(305:330),Ng(305:330),'k');
                        plot(Eg(331:356),Ng(331:356),'k');
                        plot(Eg(357:end),Ng(357:end),'k');
                    else
                        plot(Eg,Ng,'k'); hold on; end
                   %ELA
                    Eela = (rig.ELA(ELA_d(i,1):ELA_d(i,2),1) - minE)/40;
                    Nela = (rig.ELA(ELA_d(i,1):ELA_d(i,2),2) - minN)/40; Nela = (max(Ng)-Nela);
                    plot(Eela,Nela ,'k--')
                    xlim([0 G13size(1,1)]); ylim([0 G13size(1,2)])
                    set(gca,'Ydir','reverse')
                end
                
                %Plotting dots
                E = (SWE(i).utm(:,1)-min(SWE(i).utm(:,1)))/40;
                Na = SWE(i).utm(:,2)-min(SWE(i).utm(:,2)); N = (max(Na)-Na)/40;                
                if      i==1; E = E+26; N = N+85;
                    subfig = gcf;
                elseif  i==2; E = E+11; N = N+49;
                elseif  i==3; E = E+17; N = N+21; end
                
                if      strcmp(sweDOTS,'black')
                    plot(E,N,'k.', 'MarkerSize',5); hold on
                    caxis([x_min x_max]); caxis(caxis); 
                elseif  strcmp(sweDOTS,'colour')
                    scatter(E,N , 13, SWE(i).swe,'filled'); 
                    caxis([minSWE maxSWE]); caxis(caxis); 
                elseif  strcmp(sweDOTS,'sweONswe')
                    scatter(E,N , 13, SWE(i).swe,'filled'); 
                    caxis([x_min x_max]); caxis(caxis);
                elseif   strcmp(sweDOTS,'symmetric')
                    plot(E,N,'k.', 'MarkerSize',5); hold on
                        mc = max(abs([x_min x_max]));
                    caxis([-mc mc]); caxis(caxis); 
                    C = cbrewer('div', 'PRGn', 21, 'PCHIP');
                        colormap(flipud(C))
                end        
                
                %Axis Properties
                axis equal
                s.Position   = pos_axis(i,:);
                s.Box        = 'off';  axis off 
end

%Colour bar
    c = colorbar('location','eastoutside');  ylabel(c,cLabel)
    set(c,'Position',[0.88 0.2 0.03 0.55]);

%Flow direction
if ~all(isnan(topoParam.(name)(:)))
    annotation('arrow',[.15 .19],[.32 .22]) %G4
    annotation('arrow',[.39 .31],[.55 .61]) %G2
    annotation('arrow',[.81 .74],[.47 .59]) %G13
elseif all(isnan(topoParam.(name)(:)))
    annotation('arrow',[.13 .17],[.52 .42]) %G4
    annotation('arrow',[.39 .31],[.55 .61]) %G2
    annotation('arrow',[.76 .70],[.47 .59]) %G13
end  

%Glacier labels
    annotation('textbox',[.02 .02 .1 .1],'String', 'Glacier 4','EdgeColor','none')
    annotation('textbox',[.32 .02 .1 .1],'String', 'Glacier 2','EdgeColor','none')
    annotation('textbox',[.6 .02 .1 .1],'String', 'Glacier 13','EdgeColor','none')

    
% North arrow
    Narrow = imread('Narrow.jpg');
    a = axes('position',[0.855,0.82,0.12,0.12]); 
    Nshow = imshow(Narrow, Cmap);
    colormap(a,gray)    
    axis off; 

%Scale bar
    a = axes(subfig); axis off; 
    scalebar('ScaleLength', .12, 'Location',[0.9,0.95])
    annotation(gcf,'textbox',[0.8,0.83,.08,.05],...
                'String',{'2 km'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on
    annotation(gcf,'textbox',[0.725,0.83,.05,.05],...
                'String',{'0'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on

%Font size and image size           
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];

end

