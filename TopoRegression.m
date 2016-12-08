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
options.ZZ          = 2; %exclude zigzags
run MAIN

    for i = 1:3
        y       = SWE(i).swe;
        name    = char(options.glacier(i)); 
            display(['option = ',num2str(t), ', glacier = ',name]);
        X       = topo_sampled.(name);

        [MLR(t).(name), residualsMLR(t).(name)] = MLRcalval(y, X);
        MLR(t).(name).Properties.VariableNames = strcat(MLR(t).(name).Properties.VariableNames, num2str(t));
    end
end
        clear best i name X y t glacier
        
%% Export all values
G4_mlrDensity = []; G2_mlrDensity = []; G13_mlrDensity = [];
for i = 2:9
G4_mlrDensity  = [G4_mlrDensity, MLR(i).G4(:,3)];
G2_mlrDensity  = [G2_mlrDensity, MLR(i).G2(:,3)];
G13_mlrDensity = [G13_mlrDensity, MLR(i).G13(:,3)];
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
for r = 2:9
option = r;

figure(1)
for i = 1:3
    name	= char(options.glacier(i));
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

%% Plots - MLR fit for all SWE options

R2value = [];

figure(1)
for i = 1:3

    for r = 2:9
    option = r;
    name	= char(options.glacier(i));
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
        [f.(name), g.(name)] = fit(y, X,'poly1');
        p = plot(f.(name)); hold on
        set(p,'Color',options.RGB(i,:)); set(p, 'LineWidth',1.5);     
        xlabel('Measured Winter Balance (m w.e.)'); ylabel('Modelled Winter Balance (m w.e.)');
        title(name)
                axis square;    box on        
        R2value = mean([R2value g.(name).rsquare]);
    end
        b = gca; legend(b,'off');
        dim = [b.Position(1)+0.01 b.Position(2)+.37 .3 .3];
        annotation('textbox',dim,'String', ['R^2=',num2str(round(R2value,2))],'FitBoxToText','on')
    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 12 4];
filename = 'MLRfit_allLines';
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

end

%% Plots - Residuals

figure(1)
for i = 1:3

    for r = 2:9
    option = r;    name	= char(options.glacier(i));
    
    subplot(1,3,i)
    [N, edges] = histcounts(residuals(r).(name)); hold on
    plot(N,'-o')
    %hist(residuals(r).(name)); hold on
        %b = gca; 
        %dim = [b.Position(1)+0.01 b.Position(2)+.5 .3 .3];
        %[~, chi] = chi2gof(residuals(r).(name));
        %annotation('textbox',dim,'String', ['p_{\chi} = ', num2str(chi)], 'EdgeColor','none')
    end
        legend('Option 1', 'Option 2','Option 3','Option 4','Option 5',...
            'Option 6','Option 7','Option 8')
        xlim([0 40]); ylim([0 180]);
        xlabel('Residual (m w.e.)');         title(name) 
        if i == 1; ylabel('Frequency');
        else     ylabel('');        end
end
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = 'MLRresiduals_all';
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

%% Plots - Box and whisker for coeffs from density options
clf

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

aboxplot(h,'labels',options.topoVars2, ...
    'Colormap',                 options.RGB,...
    'OutlierMarkerSize',        10,...
    'OutlierMarkerEdgeColor',   [0 0 0]); % Advanced box plot
        legend('Glacier 4','Glacier 2','Glacier 13'); % Add a legend
        ylabel('Variance Explained (%)'); hold on 

        for i = 1:7
        line([i+0.5 i+0.5],[0 100],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end


        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',22) 
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
    [pearson.(glacier) Ppearson.(glacier)] = corr(X); 
    pearson.(glacier)   = round(triu(pearson.(glacier)),2); 
    Ppearson.(glacier)  = round(triu(Ppearson.(glacier)),2);
%     matrix2latex(pearson.(glacier),'/home/glaciology1/Documents/MastersDocuments/Methods/temp.txt',...
%         'rowLabels',options.topoVars_xunit,'columnLabels',options.topoVars_xunit, 'alignment','c')
end
    
%Full topo
for i = 1:3
X = [];
glacier = char(options.glacier(i));
params = fieldnames(topo_full.(glacier));
    for p = 1:length(params)
        PP = char(params(p));
        X = [X, topo_full.(glacier).(PP)(:)];
    end
    [pearson.(glacier) Ppearson.(glacier)] = corr(X); 
    pearson.(glacier)   = round(triu(pearson.(glacier)),2); 
    Ppearson.(glacier)  = round(triu(Ppearson.(glacier)),2);
end
    
%% BMS

for t = 2:9
    run OPTIONS
    options.DensitySWE  = t;
    options.ZZ          = 2; %exclude zigzags
    run MAIN
    display(['Option ', num2str(t)]);
    
    cd BMS
    BMSinit = BMS_R(SWE, topo_sampled);
    cd ..
BMS(t,:).G4     = BMSinit.G4;
BMS(t,:).G2     = BMSinit.G2;
BMS(t,:).G13     = BMSinit.G13;
end
    clear BMSinit t

% Ploting coeffs
clf
opt = 8;
II = 1; %To include intercept =1, exclude intercept =2
heads    = fieldnames(BMS(opt).G4);  heads = heads(1:end-1);
rowNames = BMS(opt).G4.Properties.RowNames(1:end-II);
    for j = 1:length(heads)
        param = char(heads(j));
        s1 = subplot(1,3,1); title('G4')
            plot(1:9,BMS(opt).G4.(param)(1:end-II),'o','MarkerSize',10); hold on
                ylabel('BMA Coefficient')
        s2 = subplot(1,3,2); title('G2')
            plot(1:9,BMS(opt).G2.(param)(1:end-II),'o','MarkerSize',10); hold on
            	ylabel('BMA Coefficient')
        s3 = subplot(1,3,3); title('G13')
            plot(1:9,BMS(opt).G13.(param)(1:end-II),'o','MarkerSize',10); hold on
                ylabel('BMA Coefficient')
    end
          legend(heads,'Location','best')
          set([s1, s2, s3],'XTick', 1:length(rowNames), 'XTickLabel',rowNames); 
          axes(s1); rotateticklabel(s1,45); 
          axes(s2); rotateticklabel(s2,45);   
          axes(s3); rotateticklabel(s3,45);

        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 7];
filename = 'BMScoeff_compare';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')        

%% Predicting

for t = 2:9
    %check coeff order - MLR
    mlrCoeff = MLR(t).G4.Properties.RowNames(1:end-2);   topoCoeff = fieldnames(topo_full_ns.G4);
    if ~isequal(mlrCoeff, topoCoeff)
        disp('Different order of coefficients between MLR and topo'); return; end
    %check coeff order - BMS
    bmaCoeff = BMS(t).G4.Properties.RowNames(1:end-2);   topoCoeff = fieldnames(topo_full_ns.G4);
    if ~isequal(bmaCoeff, topoCoeff)
        disp('Different order of coefficients between BMS and topo'); return; end
    
    for g = 1:3
    glacier = char(options.glacier(g));
        %MLR
         %Intercept
        sweMLR(t).(glacier) = repmat(MLR(t).(glacier){9,1}, size(topo_full_ns.(glacier).centreD));
         %multiply coeffs and add them
        for n = 1:length(mlrCoeff)
            param = char(mlrCoeff(n));
            sweT = topo_full_ns.(glacier).(param)*MLR(t).(glacier){n,1};
            sweMLR(t).(glacier) = sweMLR(t).(glacier) + sweT;
        end
        
        %BMS
         %Intercept
        sweBMS(t).(glacier) = repmat(BMS(t).(glacier){9,1}, size(topo_full_ns.(glacier).centreD));
         %multiply coeffs and add them
        for n = 1:length(bmaCoeff)
            param = char(bmaCoeff(n));
            sweT = topo_full_ns.(glacier).(param)*BMS(t).(glacier){n,1};
            sweBMS(t).(glacier) = sweBMS(t).(glacier) + sweT;
        end
    end
end

