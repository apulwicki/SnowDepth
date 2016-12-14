function [ ] = PlotTopoParameter( topoParam, paramName, cLabel, SWE )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
run OPTIONS.m
clf
%get colour min max
x_min = nanmin([topoParam.G4(:);topoParam.G2(:);topoParam.G13(:)]);
x_max = nanmax([topoParam.G4(:);topoParam.G2(:);topoParam.G13(:)]);
G13size = size(topoParam.G13);

pos_axis    = [0        -0.15   0.4     1.3; 
               0.19     -0.15   0.4     1.3; 
               0.50     -0.15   0.4     1.3];

%phasemap needs
%paramName = 'elevation';
if strcmp(paramName, 'aspect') || strcmp(paramName, 'northness')
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
                h = imagesc(data); hold on
                caxis([x_min max([x_max,max(data)])]); 
                set(h,'alphadata',~isnan(data))
                
                E = (SWE(i).utm(:,1)-min(SWE(i).utm(:,1)))/40;
                Na = SWE(i).utm(:,2)-min(SWE(i).utm(:,2)); N = (max(Na)-Na)/40;                
                if      i==1; E = E+25; N = N+78;
                    subfig = gcf;
                elseif  i==2; E = E+10; N = N+47;
                elseif  i==3; E = E+17; N = N+21; 
                end
                plot(E,N,'k.'); hold on
                
                axis equal
                s.Position   = pos_axis(i,:);
                s.Box        = 'off';  axis off 
end

%Colour bar
    c = colorbar('location','eastoutside');  ylabel(c,cLabel)
    set(c,'Position',[0.9 0.13 0.03 0.6])

% North arrow
    Narrow = imread('Narrow.jpg');
    a = axes('position',[0.855,0.8,0.12,0.12]); 
    Nshow = imshow(Narrow, Cmap);
    colormap(a,gray)    
    axis off; 

%Scale bar
    a = axes(subfig); axis off; 
    scalebar('ScaleLength', .12, 'Location',[0.9,0.93])
    annotation(gcf,'textbox',[0.81,0.817,.05,.05],...
                'String',{'2 km'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on
    annotation(gcf,'textbox',[0.725,0.817,.05,.05],...
                'String',{'0'}, 'LineStyle','none','FitBoxToText','off','EdgeColor',[1 1 1],'BackgroundColor',[1 1 1]); hold on

%Font size and image size           
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];

%Save figure
    filename = ['Map_',paramName];
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
%%
end

