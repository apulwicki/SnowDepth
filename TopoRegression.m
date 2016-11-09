%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

load TopoMLR.mat
%run Import_Topo.m

%% MLR - Topo Regression

% %remove aspect
% glacier = {'G4','G2','G13'};
% for i = 1:3
%    name     = char(glacier(i));
%    topo_sampled.(name) = rmfield(topo_sampled.(name),'aspect');
% end

for t = 2:9
run OPTIONS.m
options.DensitySWE  = t;
options.ZZ          = 2; %no zigzags
run MAIN

    for i = 1:3
        glacier = [4,2,13];
        y       = SWE(i).swe;
        name    = ['G', num2str(glacier(i))]; 
            display(['option = ',num2str(t), ', glacier = ',name]);
        X       = topo_sampled.(name);

        [mlr(t).(name), residuals(t).(name)] = MLRcalval(y, X);
        mlr(t).(name).Properties.VariableNames = strcat(mlr(t).(name).Properties.VariableNames, num2str(t));
    end
end
        clear best i name X y t glacier

%% Export all values
G4_mlrDensity = []; G2_mlrDensity = []; G13_mlrDensity = [];
for i = 2:9
G4_mlrDensity  = [G4_mlrDensity, mlr(i).G4(:,3)];
G2_mlrDensity  = [G2_mlrDensity, mlr(i).G2(:,3)];
G13_mlrDensity = [G13_mlrDensity, mlr(i).G13(:,3)];
end
%     writetable(G4_mlrDensity,'/home/glaciology1/Documents/Data/GlacierTopos/G4_mlrDensity.xlsx')
%     writetable(G2_mlrDensity,'/home/glaciology1/Documents/Data/GlacierTopos/G2_mlrDensity.xlsx')
%     writetable(G13_mlrDensity,'/home/glaciology1/Documents/Data/GlacierTopos/G13_mlrDensity.xlsx')
%     
    writetable(G4_mlrDensity,'/Users/Alexandra/Documents/SFU/Data/G4_mlrDensity.xlsx')
    writetable(G2_mlrDensity,'/Users/Alexandra/Documents/SFU/Data/G2_mlrDensity.xlsx')
    writetable(G13_mlrDensity,'/Users/Alexandra/Documents/SFU/Data/G13_mlrDensity.xlsx')
    
    
%Check normality of reisduals
%plotResiduals(lm.G4) %many options for plotting

clear A box centreline d distance f* G* h i min_dist NoZZ param* RGB X1 Y1
%% Plots - MLR

% Actual vs fitted data  
glacier = [4,2,13];

option = 2;

figure
for i = 1:3
    name    = ['G', num2str(glacier(i))];
    y       = SWE(i).swe;
    X       = topo_sampled.(name);
        params = fieldnames(X); 
        M = X.(char(params(1)));
        for j = 2:length(params)
            field = char(params(j));
            M = [M, X.(field)];
        end
fitted_swe = sum(repmat(mlr(option).(name){:,1}',length(M),1).*[ones(length(M),1),M],2);

    subplot(1,3,i)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
        plot(y, fitted_swe, '.', 'Color', options.RGB(i,:),'MarkerSize',13); hold on

        [f.(name), g.(name)] = fit(y, fitted_swe,'poly1');
        p = plot(f.(name)); hold on
        set(p,'Color',RGB(i,:)); set(p, 'LineWidth',1.5);     
        xlabel('Original SWE (m)'); ylabel('MLR SWE (m)');
        legend('Reference Line',['Glacier ',num2str(glacier(i))],['R^2=',num2str(round(g.(name).rsquare,2))])
   
end

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = 'MLRfit';
% print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

%% Plots - Box and whisker for coeffs from density options

%Rearrange to compare density options
params = mlr(2).G4.Properties.RowNames;
box.G4 = []; box.G2 = []; box.G13 = [];
for i = 2:9
    box.G4 = [box.G4; table2array(mlr(i).G4(:,3))'];
    box.G2 = [box.G2; table2array(mlr(i).G2(:,3))'];
    box.G13 = [box.G13; table2array(mlr(i).G13(:,3))'];
end
%No intercept option
box.G4 = box.G4(:,2:end); box.G2 = box.G2(:,2:end); box.G13 = box.G13(:,2:end); params = params(2:end);

h = cat(1, reshape(box.G4,[1 size(box.G4)]), reshape(box.G2,[1 size(box.G2)]),...
            reshape(box.G13,[1 size(box.G13)]));

    RGB = [9, 132, 103; 224, 187, 2; 130, 75, 135]/255;
aboxplot(h,'labels',params, ...
    'Colormap',                 RGB,...
    'OutlierMarkerSize',        10,...
    'OutlierMarkerEdgeColor',   [0 0 0]); % Advanced box plot
        legend('Glacier 4','Glacier 2','Glacier 13'); % Add a legend
        ylabel('% Variance Explained')
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 10];

filename = 'Coeffs_DensityOpts';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')

%% ANOVA between topographic params

for i = 1:3
    y  = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    topo    = [aspect.(name); northness.(name); profileCurve.(name); ...
                tangentCurve.(name); slope.(name); elevation.(name); Sx.(name)];
    group   = cellstr([repmat('aspect',length(aspect.(name)),1); repmat('northn',length(northness.(name)),1);...
                repmat('profCu',length(profileCurve.(name)),1); repmat('tangCu',length(tangentCurve.(name)),1);...
                repmat('slopee',length(slope.(name)),1); repmat('elevat',length(elevation.(name)),1); repmat('Sxxxxx',length(Sx.(name)),1)]);
    [p,tbl,stats] = anova1(topo,group);
        
            
end    
    
%% Range of params sampled

%run Import_Topo.m
    
header = fieldnames(topo_full.G4);
units   = {'', '(^{\circ})','(m a.s.l)','','(m^{-1})','(^{\circ})','(m^{-1})','(m)'};
glacier = {'G4','G2','G13'}; N = zeros(3,10); edges = zeros(3,11); Nall = N; edgesall = edges;
for r = 1:length(header)
figure
    param = char(header(r));
    for i = 1:3
        name    = char(glacier(i)); 
        [N(i,:), edges(i,:)] = histcounts(topo_sampled_ns.(name).(param),10); 
        [Nall(i,:), edgesall(i,:)] = histcounts(topo_full.(name).(param),10); 
    end
    for i = 1:3
        name = char(glacier(i)); 
        a(i) = subplot(1,3,i);
            plot((edges(i,1:end-1)+edges(i,2:end))/2,N(i,:)); hold on %sampled values
            plot((edgesall(i,1:end-1)+edgesall(i,2:end))/2,Nall(i,:))
            xlabel({[char(header(r)), ' ', char(units(r))], name});     ylabel('Freq.')
            axis tight
            if i == 1; legend('Sampled','Full range'); end
    end
    linkaxes(flip(a));
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = ['SampledRangeTopo_',header{r}];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
end 

    close all
    clear v i r name header glacier N* a filename edges* param fig units


