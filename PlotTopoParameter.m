function [ ] = PlotTopoParameter( topoParam, paramName, cLabel, SWE, sweDOTS, massB)
%Inputs: data structure, parameter name, colorbar label, full SWE, color of
%dots, 'massB'
%Color of dots: black, colour, sweONswe, symmetric
global options
    rig = options.rig;
    ELA_d = [1 7; 16 23; 8 15];
clf
colormap('default')

%Set up SWE structure if not corretly formated
if length(SWE) == 1
   SWEtemp(1).swe = SWE.G4(:,1);    SWEtemp(2).swe = SWE.G2(:,1);    SWEtemp(3).swe = SWE.G13(:,1); 
   SWEtemp(1).utm = SWE.G4(:,2:3);  SWEtemp(2).utm = SWE.G2(:,2:3);  SWEtemp(3).utm = SWE.G13(:,2:3); 
   SWE = SWEtemp;
end

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
    x_max = 60;    maxSWE  = x_max;
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

%other color schemes
if  strcmp(paramName, 'uncertainity')
    Cmap = cbrewer('seq','OrRd',50,'PCHIP');
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
                
                    minE = min(rig.(name)(:,1));
                        minN = min(rig.(name)(:,2));
                    Eg = (rig.(name)(:,1) - minE)/40;
                    Ng = (rig.(name)(:,2) - minN)/40;  Ng = (max(Ng)-Ng);
                    Ng = Ng + (G13size(1,1)-max(Ng));
                    
                if all(isnan(topoParam.(name)(:))) || strcmp(sweDOTS,'symmetric')
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
                 E = (SWE(i).utm(:,1)-minE)/40;
                 Na = (SWE(i).utm(:,2)-minN)/40; N = max(Ng)-Na;
                if      strcmp(sweDOTS,'black')
                    plot(E,N,'k.', 'MarkerSize',5); hold on
                    caxis([x_min x_max]); caxis(caxis); 
                elseif  strcmp(sweDOTS,'colour')
                    INN = ~isnan(SWE(i).swe);
                    scatter(E(INN),N(INN), 13, SWE(i).swe(INN),'filled'); 
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

% Glacier labels
    annotation('textbox',[.02 .05 .1 .1],'String', 'Glacier 4','EdgeColor','none','FontWeight','bold')
    annotation('textbox',[.30 .05 .1 .1],'String', 'Glacier 2','EdgeColor','none','FontWeight','bold')
    annotation('textbox',[.59 .05 .1 .1],'String', 'Glacier 13','EdgeColor','none','FontWeight','bold')

% Winter balance
     if strcmp(massB, 'massB')
     annotation('textbox',[.02 0 .1 .1],'String',...
         [num2str(round(nanmean(topoParam.G4(:)),2), '%.2f'),' m w.e.'],'EdgeColor','none')    
     annotation('textbox',[.30 0 .1 .1],'String',...
         [num2str(round(nanmean(topoParam.G2(:)),2), '%.2f'),' m w.e.'],'EdgeColor','none')    
     annotation('textbox',[.59 0 .1 .1],'String',...
         [num2str(round(nanmean(topoParam.G13(:)),2), '%.2f'),' m w.e.'],'EdgeColor','none')    
     end
% North arrow
    Narrow = imread('Narrow.jpg');
    a = axes('position',[0.855,0.82,0.12,0.12]); 
    imshow(Narrow, Cmap);
    colormap(a,gray)    
    axis off; 

% Scale bar
    %a = axes(subfig); 
    axes('position',[0.7,0.82,0.2,0.12]); axis off; 
    scalebar('ScaleLength', 0.4, 'Location',[0.58,0.55])
    annotation(gcf,'textbox',[0.8,0.83,.08,.05],...
                'String',{'2 km'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on
    annotation(gcf,'textbox',[0.727,0.83,.05,.05],...
                'String',{'0'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on

% Font size and image size           
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];

end

