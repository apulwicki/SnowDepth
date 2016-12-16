%% %%%%%%%%%%%%%%%%%%%%% BMS and MLR Plots %%%%%%%%%%%%%%%%%%%%%

% data            = BMS;
% residualsdata   = residualsBMS;
% type            = 'BMS';

data            = MLR;
residualsdata   = residualsMLR;
type            = 'MLR';

%% Plotting - fit to observed

% Actual vs fitted data  
% for r = 2:9
option = 8;
r = 8;

figure(1)
for i = 1:3
    glacier     = char(options.glacier(i));
    yObserved   = SWE(i).swe;
    X           = struct2table(topo_sampled.(glacier));	X = X{:,:};
    coeffs      = repmat(data(option).(glacier){1:end-2,1}',length(X),1);   
    yModel      = sum(X.*coeffs,2) + data(option).(glacier){end-1,1};      
  
    subplot(1,3,i)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
        plot(yObserved, yModel, '.', 'Color', options.RGB(i,:),'MarkerSize',13); hold on

        [f.(glacier), g.(glacier)] = fit(yObserved, yModel,'poly1');
        p = plot(f.(glacier)); hold on
        set(p,'Color',options.RGB(i,:)); set(p, 'LineWidth',1.5);     
        xlabel('Measured Winter Balance (m w.e.)'); ylabel('Modelled Winter Balance (m w.e.)');
        title(['Glacier ',glacier(2:end)])
                axis square;    box on
        b = gca; legend(b,'off');
        dim = [b.Position(1)+0.01 b.Position(2)+.37 .3 .3];
        annotation('textbox',dim,'String', ['R^2=',num2str(round(g.(glacier).rsquare,2))],'FitBoxToText','on')
end

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 12 4];
filename = [type,'fit_opt',num2str(r)];
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

% clf
% end
    clear b coeffs dim f fig filename g glacier i line option p r X yModel yObserved
%% Plots - BMS Box and whisker for coeffs from density options
clf

%Rearrange to compare density options
params = data(2).G4.Properties.RowNames(1:end-2);
box.G4 = []; box.G2 = []; box.G13 = [];
for i = 2:9
    box.G4  = [box.G4;  table2array(data(i).G4(1:end-2,2))'];
    box.G2  = [box.G2;  table2array(data(i).G2(1:end-2,2))'];
    box.G13 = [box.G13; table2array(data(i).G13(1:end-2,2))'];
end

% box.G4 = log(box.G4); box.G2 = log(box.G2); box.G13 = log(box.G13);

h = cat(1, reshape(box.G4,[1 size(box.G4)]), reshape(box.G2,[1 size(box.G2)]),...
            reshape(box.G13,[1 size(box.G13)]));

aboxplot(h,'labels',options.topoVars, ...
    'Colormap',                 options.RGB,...
    'OutlierMarkerSize',        10,...
    'OutlierMarkerEdgeColor',   [0 0 0]); % Advanced box plot
        legend('Glacier 4','Glacier 2','Glacier 13'); % Add a legend
        ylabel('Variance Explained (%)'); hold on 

        for i = 1:7
        line([i+0.5 i+0.5],[0 70],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end


        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',22) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 10];

filename = [type,'Coeffs_DensityOpts'];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
    clear box fig filename h i params
%% Plots - BMS fit for all SWE options
clf
R2value = [];

figure(1)
for i = 1:3

    for r = 2:9
    option = r;
    name	= char(options.glacier(i));
    y       = SWE(i).swe;
    X       = topo_sampled.(name);
        params = data(option).(name).Properties.RowNames(1:end-2); 
        for j = 1:length(params)
            field       = char(params(j));
            X.(field)   = X.(field)*data(option).(name){j,1};
        end
	X       = struct2table(X);
    X       = sum(X{:,:},2)+data(option).(name){end-1,1};

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
filename = [type,'fit_allLines'];
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

end
        clear b p R2value dim f g j r fig filename h i params X y

%% Plots - Residuals
clf
for i = 1:3

    for r = 2:9
    option = r;    name	= char(options.glacier(i));
    
    subplot(1,3,i)
    [N, ~] = histcounts(residualsdata(r).(name)); hold on
    plot(N,'-o')
    %hist(residuals(r).(name)); hold on
        %b = gca; 
        %dim = [b.Position(1)+0.01 b.Position(2)+.5 .3 .3];
        %[~, chi] = chi2gof(residuals(r).(name));
        %annotation('textbox',dim,'String', ['p_{\chi} = ', num2str(chi)], 'EdgeColor','none')
    end
        legend('Option 1', 'Option 2','Option 3','Option 4','Option 5',...
            'Option 6','Option 7','Option 8')
        xlim([0 40]); ylim([0 190]);
        xlabel('Residual (m w.e.)');         title(name) 
        if i == 1; ylabel('Frequency');
        else     ylabel('');        end
end
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
filename = [type,'residuals_all'];
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')
    clear edges field line N name option
    
%% %%%%%%%%%%%%%%%%%%%%% OTHER %%%%%%%%%%%%%%%%%%%%%

%% Distirbution of # of params from all possible combos 
n = 8;
C = logical(dec2bin(0:(2^n)-1)=='1');       C = C(2:end,:);
C = double(C);                              C = sum(C,2);

hist(C,n)
    xlabel('Number of parameters')
    ylabel('Number of models')

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 6 6];
filename = 'DistributionOfNumParams_topoRegress';
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')

    clear C n fig filename
    
%% Range of params sampled

%run Import_Topo.m
    
header = fieldnames(topo_full.G4);
for r = 1:length(header)
figure
    param = char(header(r));
    for i = 1:3
        name = char(options.glacier(i)); 
        a(i) = subplot(1,3,i);
            h(i).f = histogram(topo_full_ns.(name).(param)(:),25); hold on
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
    clear v i r name header glacier N* a filename edges* param fig units dim h 
    
%% Ploting - MLR and BMS coeffs
clf
II = 2; %Include intercept =1, exlude =2
%Rearrange to compare density options
params = BMS(2).G4.Properties.RowNames(1:end-II);
box.G4 = []; box.G2 = []; box.G13 = [];
for i = 2:9
    box.G4  = [box.G4;  table2array(BMS(i).G4(1:end-II,1))'; table2array(MLR(i).G4(1:end-II,1))'];
    box.G2  = [box.G2;  table2array(BMS(i).G2(1:end-II,1))'; table2array(MLR(i).G2(1:end-II,1))'];
    box.G13 = [box.G13; table2array(BMS(i).G13(1:end-II,1))'; table2array(MLR(i).G13(1:end-II,1))'];
end
h = cat(1, reshape(box.G4,[1 size(box.G4)]), reshape(box.G2,[1 size(box.G2)]),...
            reshape(box.G13,[1 size(box.G13)]));

aboxplot(h,'labels',options.topoVars, ...
    'Colormap',                 options.RGB,...
    'OutlierMarkerSize',        10,...
    'OutlierMarkerEdgeColor',   [0 0 0]); % Advanced box plot
        legend('Glacier 4','Glacier 2','Glacier 13'); % Add a legend
        ylabel('Coefficient Value'); hold on 
    
        for i = 1:7
        line([i+0.5 i+0.5],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        line(xlim,[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end
        
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',22) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 10];
filename = 'Coeff_compare';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')        
     


    clear box fig filename h i II params