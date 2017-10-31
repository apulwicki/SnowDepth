function [ ] = PlotTopoParameter_DD( topoParam, cLabel, P1, P2, sweDOTS )
%Used for plotting sampling pattern or overlaying multiple set of
%measurement locations on the glaciers
%Inputs: data structure, parameter name, colorbar label, pattern1, pattern2, color of
%dots, 'massB'
%Color of dots: black, colour, sweONswe, symmetric
global options
    rig = options.rig;
    ELA_d = [1 7; 16 23; 8 15];
clf
colormap('default')

%Set up SWE structure if not corretly formated
if length(P1) == 1
   SWEtemp(1).swe = P1.G4(:,1);    SWEtemp(2).swe = P1.G2(:,1);    SWEtemp(3).swe = P1.G13(:,1); 
   SWEtemp(1).utm = P1.G4(:,2:3);  SWEtemp(2).utm = P1.G2(:,2:3);  SWEtemp(3).utm = P1.G13(:,2:3); 
   P1 = SWEtemp;
   
   SWEtemp(1).swe = P2.G4(:,1);    SWEtemp(2).swe = P2.G2(:,1);    SWEtemp(3).swe = P2.G13(:,1); 
   SWEtemp(1).utm = P2.G4(:,2:3);  SWEtemp(2).utm = P2.G2(:,2:3);  SWEtemp(3).utm = P2.G13(:,2:3); 
   P2 = SWEtemp;
end

%get colour min max
     %Manual colour range = 1
   fixed_color = 1;
if fixed_color == 0 
    x_min   = nanmin([topoParam.G4(:);topoParam.G2(:);topoParam.G13(:)]);
    x_max   = nanmax([topoParam.G4(:);topoParam.G2(:);topoParam.G13(:)]);
    minSWE  = nanmin([P1(1).swe(:);P1(2).swe(:);P1(3).swe(:)]);
    maxSWE  = nanmax([P1(1).swe(:);P1(2).swe(:);P1(3).swe(:)]);
elseif fixed_color == 1
    x_min = 0;      minSWE  = x_min;
    x_max = 1;    maxSWE  = x_max;
end

G13size = size(topoParam.G13);

width  = 0.59;
height = 2.3;
ystart = -0.7;
pos_axis    = [0.00     ystart   width     height; 
               0.095     ystart   width     height; 
               0.50     ystart   width     height];

%Colormap choice
    Cmap = parula;
    %Cmap = cbrewer('div','YlGnBu',100);
 
for i = 1:3
        name    = char(options.glacier(i)); 
        Gsize   = size(topoParam.(name));
        data    = [ nan(G13size(1,1)-Gsize(1,1),Gsize(1,2)); topoParam.(name)];
        data    = [data, nan(G13size(1,1), G13size(1,2)-Gsize(1,2))];
        
        s(i) = subplot(1,3,i);
                h(i) = imagesc(data); hold on
                colormap(Cmap); set(h(i),'alphadata',~isnan(data));
                        minE = min(rig.(name)(:,1));
                        minN = min(rig.(name)(:,2));
                    Ng = (rig.(name)(:,2) - minN)/40;  Ng = (max(Ng)-Ng);
                    Ng = Ng + (G13size(1,1)-max(Ng));
                    

                   %ELA
                    Eela = (rig.ELA(ELA_d(i,1):ELA_d(i,2),1) - minE)/40;
                    Nela = (rig.ELA(ELA_d(i,1):ELA_d(i,2),2) - minN)/40; Nela = (max(Ng)-Nela);
                    plot(Eela,Nela ,'k--')
                    xlim([0 G13size(1,1)]); ylim([0 G13size(1,2)])
                    set(gca,'Ydir','reverse')
                
                %Plotting dots
                 caxis([x_min x_max]); caxis(caxis);
                 E1 = (P1(i).utm(:,1)-minE)/40;
                 Na = (P1(i).utm(:,2)-minN)/40; N1 = max(Ng)-Na;
                 E2 = (P2(i).utm(:,1)-minE)/40;
                 Na = (P2(i).utm(:,2)-minN)/40; N2 = max(Ng)-Na;
                if      strcmp(sweDOTS,'black')
                    plot(E1,N1,'k.', 'MarkerSize',1); hold on
                    plot(E2,N2,'ko', 'MarkerSize',2); 
                end        
                
                %Axis Properties
                axis equal
                ylim([0 160]) 
                s(i).Box        = 'off';  axis off 

end

for i = 1:3;
        s(i).Position   = pos_axis(i,:);
end

%Colour bar
    c = colorbar;  ylabel(c,cLabel)
    set(c,'Position',[0.87 0.5 0.03 0.37]);

%Flow direction
    arG4 = annotation('arrow',[.18 .25],[.32 .26]); arG4.HeadLength = 5; arG4.HeadWidth = 5; %G4
    arG2 = annotation('arrow',[.39 .31],[.49 .52]); arG2.HeadLength = 5; arG2.HeadWidth = 5; %G2
    arG3 = annotation('arrow',[.78 .70],[.57 .64]); arG3.HeadLength = 5; arG3.HeadWidth = 5; %G13


% Glacier labels
    annotation('textbox',[.06 0.05 .1 .1],'String', 'Glacier 4','EdgeColor','none','FontWeight','bold')
    annotation('textbox',[.40 0.05 .1 .1],'String', 'Glacier 2','EdgeColor','none','FontWeight','bold')
    annotation('textbox',[.67 0.05 .1 .1],'String', 'Glacier 13','EdgeColor','none','FontWeight','bold')

% North arrow
    Narrow = imread('Narrow.jpg');
    a = axes('position',[0,0.75,0.1,0.1]); 
    imshow(Narrow, Cmap);
    colormap(a,gray)    
    axis off; 

% Scale bar
    axes('position',[0.16,0.75,0.2,0.12]); axis off; 
    scalebar('ScaleLength', 0.8, 'Location',[0.18,0.55])
    annotation(gcf,'textbox',[0.3,0.76,0.15,0.05],...
                'String',{'2 km'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on
    annotation(gcf,'textbox',[0.17,0.76,.05,.05],...
                'String',{'0'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on

end

