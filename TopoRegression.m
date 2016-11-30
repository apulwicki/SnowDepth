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
options.ZZ          = 1; %include zigzags
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

% Sort topo params

order = mlr(2).G2.Properties.RowNames(2:end);
for i = 1:3
   name     = char(options.glacier(i));
   topo_full.(name)         = orderfields(topo_full.(name),order);
   topo_sampled.(name)      = orderfields(topo_sampled.(name),order);
   topo_sampled_ns.(name)   = orderfields(topo_sampled_ns.(name),order);
end
        
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

for r = 2:9
option = r;

figure(1)
for i = 1:3
    name	= ['G', num2str(glacier(i))];
    y       = SWE(i).swe;
    X       = topo_sampled.(name);
        params = mlr(option).(name).Properties.RowNames; 
        for j = 2:length(params)
            field       = char(params(j));
            X.(field)   = X.(field)*mlr(option).(name){j,1};
        end
	X       = struct2table(X);
    X       = sum(X{:,:},2)+mlr(option).(name){1,1};

    subplot(1,3,i)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
        plot(y, X, '.', 'Color', options.RGB(i,:),'MarkerSize',13); hold on

        [f.(name), g.(name)] = fit(y, X,'poly1');
        p = plot(f.(name)); hold on
        set(p,'Color',options.RGB(i,:)); set(p, 'LineWidth',1.5);     
        xlabel('Measured Winter Balance (m w.e.)'); ylabel('Modelled Winter Balance (m w.e.)');
        title(['Glacier ',num2str(glacier(i))])
                axis square;    box on
        b = gca; legend(b,'off');
        dim = [b.Position(1)+0.01 b.Position(2)+.37 .3 .3];
        annotation('textbox',dim,'String', ['R^2=',num2str(round(g.(name).rsquare,2))],'FitBoxToText','on')
end

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 12 4];
filename = ['MLRfit_opt',num2str(r)];
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

clf
end

%% Plots - Residuals

for r = 2:9
option = r;

figure(1)
for i = 1:3
    name	= char(options.glacier(i));
    
    subplot(1,3,i)
    hist(residuals(r).(name))
        xlabel('Residual (m w.e.)'); ylabel('Frequency');
        title(name) 
        b = gca; 
        dim = [b.Position(1)+0.01 b.Position(2)+.5 .3 .3];
        [~, chi] = chi2gof(residuals(r).(name));
        annotation('textbox',dim,'String', ['p_{\chi} = ', num2str(chi)], 'EdgeColor','none')
end

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = ['MLRresiduals_opt',num2str(r)];
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

clf
end


%% Plots - Box and whisker for coeffs from density options

%Rearrange to compare density options
params = mlr(2).G4.Properties.RowNames;
box.G4 = []; box.G2 = []; box.G13 = [];
for i = 2:9
    box.G4 = [box.G4; table2array(mlr(i).G4(:,2))'];
    box.G2 = [box.G2; table2array(mlr(i).G2(:,2))'];
    box.G13 = [box.G13; table2array(mlr(i).G13(:,2))'];
end

% box.G4 = log(box.G4); box.G2 = log(box.G2); box.G13 = log(box.G13);

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
        ylabel('% Variance Explained'); hold on 

%         plot(ylim,[1 1])

%         hax=axes; SP=1; %your point goes here 
%         line([SP SP],get(hax,'YLim'),'Color',[1 0 0])

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
for r = 1:length(header)
figure
    param = char(header(r));
    for i = 1:3
        name = char(options.glacier(i)); 
        a(i) = subplot(1,3,i);
            h(i).f = histogram(topo_full.(name).(param)(:),25); hold on
            h(i).s = histogram(topo_sampled_ns.(name).(param),25); 
            
            dim = [a(i).Position(1)+0.01 a(i).Position(2)*5.41 .3 .3];
            annotation('textbox',dim,'String', name,'EdgeColor','none','FontWeight','bold')
            xlabel(char(options.topoVars(r)));     ylabel('Frequency')
            axis tight
            legend('Full range','Sampled');
    end
    set([h(1).f h(1).s h(2).f h(2).s h(3).s], 'BinEdges', h(3).f.BinEdges)
    linkaxes(flip(a));
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = ['SampledRangeTopo_',header{r}];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
end 

    close all
    clear v i r name header glacier N* a filename edges* param fig units


%% Correlation between topographic parameters 

%Sampled
for i = 1:3
X = [];
glacier = char(options.glacier(i));
params = fieldnames(topo_sampled.(glacier));
    for p = 1:length(params)
        PP = char(params(p));
        X = [X, topo_sampled.(glacier).(PP)];
    end
    pearson.(glacier) = corr(X); pearson.(glacier) = round(triu(pearson.(glacier)),2);
    matrix2latex(pearson.(glacier),'/home/glaciology1/Documents/MastersDocuments/Methods/temp.txt',...
        'rowLabels',options.topoVars_xunit,'columnLabels',options.topoVars_xunit, 'alignment','c')
end
    
    
    
    