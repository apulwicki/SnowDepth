%% %%%%%%%%%%%%%%%%%%%%% BMS and MLR Plots %%%%%%%%%%%%%%%%%%%%%
% % % 
data            = BMS;
residualsdata   = residualsBMS;
type            = 'BMS';
% 
% data            = MLR;
% residualsdata   = residualsMLR;
% type            = 'MLR';

%% Plotting - fit to observed

% Actual vs fitted data  
%  for r = 2:9
r = 8;
 option = r;

figure(1)
for i = 1:3
    glacier     = char(options.glacier(i));
    yObserved   = SWE(i).swe;
    X           = struct2table(topo_sampled.(glacier));	X = X{:,:};
    coeffs      = repmat(data(option).(glacier){1:end-3,1}',length(X),1);   
    yModel      = sum(X.*coeffs,2) + data(option).(glacier){end-2,1};      
  
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
saveFIG([type,'fit_opt',num2str(r)])

 clf
%  end
    clear b coeffs dim f fig filename g glacier i line option p r X yModel yObserved
%% Plots - Box and whisker for coeffs/semi partial from density options
clf

val = 1; Lval = 'coeff';  %coeff value
%  val = 2; Lval = 'semiR2'; %semi-partial correlation

%Rearrange to compare density options
params = data(2).G4.Properties.RowNames(1:end-3);
box.G4 = []; box.G2 = []; box.G13 = [];
for i = 2:9
    box.G4  = [box.G4;  table2array(data(i).G4(1:end-3,val))'];
    box.G2  = [box.G2;  table2array(data(i).G2(1:end-3,val))'];
    box.G13 = [box.G13; table2array(data(i).G13(1:end-3,val))'];
end

% box.G4 = log(box.G4); box.G2 = log(box.G2); box.G13 = log(box.G13);

h = cat(1, reshape(box.G4,[1 size(box.G4)]), reshape(box.G2,[1 size(box.G2)]),...
            reshape(box.G13,[1 size(box.G13)]));

aboxplot(h,'labels',options.topoVars, ...
    'Colormap',                 options.RGB,...
    'OutlierMarkerSize',        10,...
    'OutlierMarkerEdgeColor',   [0 0 0]); % Advanced box plot
        legend('Glacier 4','Glacier 2','Glacier 13'); % Add a legend
        
        if val==1
            ylabel('Regression coefficient'); hold on 
            line(xlim,[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        elseif val==2
            ylabel('Semi-partial R^2'); hold on; end 

        for i = 1:7
        line([i+0.5 i+0.5],ylim,'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end


        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',22) 
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 10];

saveFIG([type,Lval,'_DensityOpts'])
     clear box fig filename h i params
%% Plots - fit for all SWE options
clf
R2value = [];

figure(1)
for i = 1:3

    for r = 2:9
    option = r;
    name	= char(options.glacier(i));
    y       = SWE(i).swe;
    X       = topo_sampled.(name);
        params = data(option).(name).Properties.RowNames(1:end-3); 
        for j = 1:length(params)
            field       = char(params(j));
            X.(field)   = X.(field)*data(option).(name){j,1};
        end
	X       = struct2table(X);
    X       = sum(X{:,:},2)+data(option).(name){end-2,1};

    subplot(1,3,i)
    axis([0 1.2 0 1.2]);    line = refline(1,0);    line.Color = 'k'; line.LineStyle = '--'; hold on
        [f.(name), g.(name)] = fit(y, X,'poly1');
        p = plot(f.(name)); hold on
        set(p,'Color',options.RGB(i,:)); set(p, 'LineWidth',1.5);     
        xlabel('Measured Winter Balance (m w.e.)'); ylabel('Modelled Winter Balance (m w.e.)');
        title(name)
                axis square;    box on        
        R2value = mean([R2value g.(name).rsquare]);
        if r == 8; p.Color = [0 0 0]; end
    end
        b = gca; legend(b,'off');
        dim = [b.Position(1)+0.01 b.Position(2)+.37 .3 .3];
        annotation('textbox',dim,'String', ['R^2=',num2str(round(R2value,2))],'FitBoxToText','on')
    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 12 4];
saveFIG([type,'fit_allLines'])

end
        clear b p R2value dim f g j r fig filename h i params X y

%% Plots - Residuals
clf
for i = 1:3

    for r = 2:9
    option = r;    name	= char(options.glacier(i));
    
    subplot(1,3,i)
    [N, edges] = histcounts(residualsdata(r).(name)); hold on
    edges = mean([edges(1:end-1); edges(2:end)],1);
    plot(edges,N,'-', 'LineWidth',2)
    %ax = gca; ax.XTick = edges;
    %hist(residuals(r).(name)); hold on
        %b = gca; 
        %dim = [b.Position(1)+0.01 b.Position(2)+.5 .3 .3];
        %[~, chi] = chi2gof(residuals(r).(name));
        %annotation('textbox',dim,'String', ['p_{\chi} = ', num2str(chi)], 'EdgeColor','none')
    end
        legend('Option 1', 'Option 2','Option 3','Option 4','Option 5',...
            'Option 6','Option 7','Option 8')
        xlim([-0.6 0.6]); ylim([0 180]);
        xlabel('Residual (m w.e.)');         title(name) 
        if i == 1; ylabel('Frequency');
        else     ylabel('');        end
end
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 15 5];
filename = [type,'residuals_all'];
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')
    clear edges field line N name option fig filename i r 
    
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
            h(i).f = histogram(topo_full_ns.(name).(param)(:),25,'FaceColor','w'); hold on
            h(i).s = histogram(topo_sampled_ns.(name).(param),25,'FaceColor','k'); 
            
            hist_stats_full.(name)(r,:) = [nanmean(topo_full_ns.(name).(param)(:)),...
                                            nanstd(topo_full_ns.(name).(param)(:)),...
                                            skewness(topo_full_ns.(name).(param)(:)),...
                                            kurtosis(topo_full_ns.(name).(param)(:))];
            hist_stats_sample.(name)(r,:) = [nanmean(topo_sampled_ns.(name).(param)(:)),...
                                            nanstd(topo_sampled_ns.(name).(param)(:)),...
                                            skewness(topo_sampled_ns.(name).(param)(:)),...
                                            kurtosis(topo_sampled_ns.(name).(param)(:))];                            
                                        
            dim = [a(i).Position(1)+0.01 a(i).Position(2)*5.41 .3 .3];
            annotation('textbox',dim,'String', name,'EdgeColor','none','FontWeight','bold')
            xlabel(char(options.topoVarsUnits(r)));     ylabel('Frequency')
            axis tight
            legend('Full range','Sampled');
    end
    set([h(1).f h(1).s h(2).f h(2).s h(3).s], 'BinEdges', h(3).f.BinEdges)
    linkaxes(flip(a));
        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
% saveFIG(['SampledRangeTopo_',header{r}])
end 

    close all
    clear v i r name header glacier N* a filename edges* param fig units dim h 
    
%% Plotting - MLR and BMS coeffs
II = 3; %Include intercept =2, exlude = 3
%Rearrange to compare density options
params = BMS(2).G4.Properties.RowNames(1:end-II);
for g = 1:3
    boxtemp.MLR = []; boxtemp.BMS = [];
    glacier = char(options.glacier(g)); 
for i = 2:9
    boxtemp.MLR  = [boxtemp.MLR;  table2array(MLR(i).(glacier)(1:end-II,1))']; 
    boxtemp.BMS  = [boxtemp.BMS;  table2array(BMS(i).(glacier)(1:end-II,1))'];
    h.(glacier) = cat(1, reshape(boxtemp.MLR,[1 size(boxtemp.MLR)]),...
                         reshape(boxtemp.BMS,[1 size(boxtemp.BMS)]));
end        
end        

clf
figure(1)
fig=gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 4];
for g = 1:3
    glacier = char(options.glacier(g));
s(1) = subplot(1,3,g)
aboxplot(h.(glacier),'labels',options.topoVars, ...
    'Colormap',                 [144 195 212; 245 177 29]/255,...
    'OutlierMarkerSize',        10,...
    'OutlierMarkerEdgeColor',   [0 0 0]); % Advanced box plot
        ylabel('Coefficient Value'); hold on 
        legend('MLR','BMA', 'Location','north'); % Add a legend
        
        YY = ylim
        for i = 1:6
                ylim_temp = YY;
        line([i+0.5 i+0.5],[ylim_temp(1,1)-0.01 ylim_temp(1,2)+0.01],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        line(xlim,[0 0],'Color',[0 0 0],'LineStyle','--','LineWidth', 0.5)
        end

end

        annotation('textbox',[0.21 0.8 0.2 0.2],'String', 'Glacier 4', 'EdgeColor','none')
        annotation('textbox',[0.49 0.8 0.2 0.2],'String', 'Glacier 2', 'EdgeColor','none')
        annotation('textbox',[0.76 0.8 0.2 0.2],'String', 'Glacier 13', 'EdgeColor','none')
       set(findall(fig,'-property','FontSize'),'FontSize',15)
saveFIG('CoeffBoxplot_BMSMLRcompare')    

    clear boxtemp fig filename h i II params g glacier
    

%% Stats for predicted SWE

%  method = 'BMS';
% res = residualsBMS(8);
method = 'MLR';
res = residualsMLR(8);


for g = 1:3
    glacier = char(options.glacier(g));
     display([glacier,' ', num2str(round(nanmean(nanmean(stackSWE.(method).(glacier)(:),3)),2)),' ',...
                num2str(round(nanmin(nanmin(stackSWE.(method).(glacier)(:),[],3)),2)),' ',...
                num2str(round(nanmax(nanmax(stackSWE.(method).(glacier)(:),[],3)),2))])
end    

clf
 %Boxplot of estimated SWE
    T = [reshape(stackSWE.(method).G4(:,:,7),[],1);...
        reshape(stackSWE.(method).G2(:,:,7),[],1);...
        reshape(stackSWE.(method).G13(:,:,7),[],1)];
    G = [repmat('G04',length(reshape(stackSWE.(method).G4(:,:,7),[],1)),1); ...
         repmat('G02',length(reshape(stackSWE.(method).G2(:,:,7),[],1)),1);...
         repmat('G13',length(reshape(stackSWE.(method).G13(:,:,7),[],1)),1)];
boxplot(T,G,'Labels',{'Glacier 4','Glacier 2','Glacier 13'})
    ylabel([{'SWE estimated with'},{[method,' coefficients (m w.e)']}])
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
    saveFIG(['ModelledSWE_box_',method])
   
 %Boxplot of residuals
    T = [res.G4(:); res.G2(:); res.G13(:)];
    G = [repmat('G04',length(res.G4(:)),1); ...
         repmat('G02',length(res.G2(:)),1);...
         repmat('G13',length(res.G13(:)),1)];
boxplot(T,G,'Labels',{'Glacier 4','Glacier 2','Glacier 13'})
    ylabel([method,' Residuals (m w.e.)'])
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
     saveFIG(['residuals_box_',method])
 
%Map of residuals at sampling locations

param = 'empty';
for g = 1:3
    glacier = char(options.glacier(g));
topoParam.(glacier)  = NaN(size(topo_full_ns.(glacier).elevation));
topoParam.rig = rig;

resZ(g).swe = res.(glacier); 
resZ(g).utm = SWE(g).utm;
end

figure(3)
PlotTopoParameter(topoParam,param, 'BMS Residuals (m w.e.)', resZ, 'colour')
    C = cbrewer('div', 'RdYlBu', 20, 'PCHIP');
    colormap(flipud(C))
        saveFIG(['residualsMap_',method])