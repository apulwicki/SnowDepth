%% Point Scale

% Box plots for each glacier

 % With zigzags
run OPTIONS.m
options.ZZ = 1;
run MAIN

boxplot([SWE(1).depth;SWE(2).depth;SWE(3).depth], [SWE(1).glacier;SWE(2).glacier;SWE(3).glacier],...
    'GroupOrder',{'G04','G02','G13'},'labels',{'Glacier 4','Glacier 2','Glacier 13'})
    ylabel('Snow depth (cm)')
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',17) 
    
    print([options.path1,'box_depth_wZZ'],'-dpng'); print([options.path2,'box_depth_wZZ'],'-dpng') 
        clear fig
        
 % Remove zigzags
run OPTIONS.m
options.ZZ = 2;
run MAIN

boxplot([SWE(1).depth;SWE(2).depth;SWE(3).depth], [SWE(1).glacier;SWE(2).glacier;SWE(3).glacier],...
    'GroupOrder',{'G04','G02','G13'},'labels',{'Glacier 4','Glacier 2','Glacier 13'})
    ylabel('Snow depth (cm)')
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',17) 
    
    print([options.path1,'box_depth_noZZ'],'-dpng'); print([options.path2,'box_depth_noZZ'],'-dpng') 
        clear fig        
%% Zigzag Plot

zig_lab = ['G04\_Z3A\_ZZ0'; 'G04\_Z2A\_ZZ0'; 'G04\_Z5B\_ZZ0';...
       'G02\_Z5C\_ZZ0'; 'G02\_Z7A\_ZZ0'; 'G02\_Z3B\_ZZ0'; ...
       'G13\_Z7C\_ZZ0'; 'G13\_Z4C\_ZZ0'; 'G13\_Z3B\_ZZ0'; 'G13\_Z5A\_ZZ0'];

for j = 2%1:3
    T1 = cellstr(char(SWE(j).label));

    for i = 5%1:size(zig_lab,1)
        T2 = find(~cellfun('isempty',strfind(T1,zig_lab(i,:))));
        if ~isempty(T2)
        data = [SWE(j).swe(T2)*100, SWE(j).utm(T2,1:2)];
            data(:,2) = max(data(:,2))-data(:,2);
            data(:,3) = max(data(:,3))-data(:,3);

        scatter(data(:,2),data(:,3),30, data(:,1),'filled')
%             title(zig_lab(i,1:8))
%             str = {strcat('mean= ', num2str(round(mean(data(:,1)),2)),'cm SWE'), ...
%                 strcat('std= ', num2str(round(std(data(:,1)),2)),'cm SWE')};
%             dim = [.13 .5 .3 .3];
%             annotation('textbox',dim,'String', str,'FitBoxToText','on')
            xlabel('Distance (m)')
            ylabel('Distance (m)')
                c = colorbar;
                c.Label.String = 'SWE (cm)';
            fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',20) 
            filename = [zig_lab(i,1:3), zig_lab(i,6:8)];
            print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
            clf

        d_ZZ = variogramAlex(data, 2, 40);
        vario = variofitAlex(d_ZZ,zig_lab(i,1:8),0);
            
%             filename = [filename, 'variogram'];
%             fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 8 9];
%             print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
%             clf
        end
     end
end

clear c d data dim filename fig i j str T* zig vario
%close all

%% Whole G plot

j = 3;
lag_G = 15; 
    d_G = variogramAlex([SWE(j).swe*100, SWE(j).utm(:,1:2)], lag_G, 'default');
    vario_G = variofitAlex(d_G,char(SWE(j).glacier(1,1)),1);
                       
 filename = 'variogram_Glacier13';
    fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 8 9];
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
        clear j lag_G d* vario* fig filename

%% SWE calculation variations

zig_lab = ['G04\_Z3A\_ZZ0'; 'G04\_Z2A\_ZZ0'; 'G04\_Z5B\_ZZ0';...
       'G02\_Z5C\_ZZ0'; 'G02\_Z7A\_ZZ0'; 'G02\_Z3B\_ZZ0'; ...
       'G13\_Z7C\_ZZ0'; 'G13\_Z4C\_ZZ0'; 'G13\_Z3B\_ZZ0'; 'G13\_Z5A\_ZZ0'];
count_zig = 1;  params_zigzag = cell(10,5);
count_G = 1;    params_G = cell(10,5); params_GU = cell(10,5); params_GL = cell(10,5);

for opts = 1:9;
    run OPTIONS.m
    options.DensitySWE = opts;
    run MAIN.m
    
for j = 1:3
    T1 = cellstr(char(SWE(j).label));

    % zigzag
    for i = 1:size(zig_lab,1)
        T2 = find(~cellfun('isempty',strfind(T1,zig_lab(i,:))));
        if ~isempty(T2)
            if opts == 1
                data = [SWE(j).depth(T2), SWE(j).utm(T2,1:2)];
            else    
                data = [SWE(j).swe(T2)*100, SWE(j).utm(T2,1:2)];
            end
            data(:,2) = max(data(:,2))-data(:,2);
            data(:,3) = max(data(:,3))-data(:,3);
        
        d_ZZ = variogramAlex(data, 2, 40);
        vario = variofitAlex(d_ZZ,zig_lab(i,1:8),0);
        params_zigzag(count_zig,:) = [num2cell(options.DensitySWE), cellstr(zig_lab(i,1:8)), ...
                                 num2cell(vario.range), num2cell(vario.nugget), num2cell(vario.sill)];
            count_zig = count_zig+1;
        end
    end
     
    % full G
    lag_G = 15; 
    d_G = variogramAlex([SWE(j).swe*100, SWE(j).utm(:,1:2)], lag_G, 'default');
    vario_G = variofitAlex(d_G,char(SWE(j).glacier(1,1)),0);
    params_G(count_G,:) = [num2cell(options.DensitySWE), cellstr(char(SWE(j).glacier(1,1))), ...
                                 num2cell(vario_G.range), num2cell(vario_G.nugget), num2cell(vario_G.sill)];
                             
    % upper G 
    T1 = num2cell(char(SWE(j).pattern)); T2 = any([strcmp(T1(:,1),'U'), strcmp(T1(:,1),'U')],2);
    d_GU = variogramAlex([SWE(j).swe(T2)*100, SWE(j).utm(T2,1:2)], lag_G, 'default');
    vario_GU = variofitAlex(d_GU,[char(SWE(j).glacier(1,1)),'UPPER'],0);
    params_GU(count_G,:) = [num2cell(options.DensitySWE), cellstr(char(SWE(j).glacier(1,1))), ...
                                 num2cell(vario_GU.range), num2cell(vario_GU.nugget), num2cell(vario_GU.sill)];
    % lower G 
    T1 = num2cell(char(SWE(j).pattern)); T2 = any([strcmp(T1(:,1),'L'), strcmp(T1(:,1),'U')],2);
    d_GL = variogramAlex([SWE(j).swe(T2)*100, SWE(j).utm(T2,1:2)], lag_G, 'default');
    vario_GL = variofitAlex(d_GL,[char(SWE(j).glacier(1,1)),'LOWER'],0);
    params_GL(count_G,:) = [num2cell(options.DensitySWE), cellstr(char(SWE(j).glacier(1,1))), ...
                                 num2cell(vario_GL.range), num2cell(vario_GL.nugget), num2cell(vario_GL.sill)];
    
            count_G = count_G+1;
end
end

values_G.range = [params_G(1:3:end,1), params_G(1:3:end,3), params_G(2:3:end,3), params_G(3:3:end,3)];
values_GL.range = [params_GL(1:3:end,1), params_GL(1:3:end,3), params_GL(2:3:end,3), params_GL(3:3:end,3)];
values_GU.range = [params_GU(1:3:end,1), params_GU(1:3:end,3), params_GU(2:3:end,3), params_GU(3:3:end,3)];

values_G.nugget = [params_G(1:3:end,1), params_G(1:3:end,4), params_G(2:3:end,4), params_G(3:3:end,4)];
values_GL.nugget = [params_GL(1:3:end,1), params_GL(1:3:end,4), params_GL(2:3:end,4), params_GL(3:3:end,4)];
values_GU.nugget = [params_GU(1:3:end,1), params_GU(1:3:end,4), params_GU(2:3:end,4), params_GU(3:3:end,4)];

values_G.sill = [params_G(1:3:end,1), params_G(1:3:end,5), params_G(2:3:end,5), params_G(3:3:end,5)];
values_GL.sill = [params_GL(1:3:end,1), params_GL(1:3:end,5), params_GL(2:3:end,5), params_GL(3:3:end,5)];
values_GU.sill = [params_GU(1:3:end,1), params_GU(1:3:end,5), params_GU(2:3:end,5), params_GU(3:3:end,5)];

values_G.all = params_G;
values_GL.all = params_GL(1:3:end,1);
values_GU.all = params_GU(1:3:end,1);


values_ZZ.range = [params_zigzag(1:10:end,1),params_zigzag(1:10:end,3),params_zigzag(2:10:end,3), ...
    params_zigzag(3:10:end,3),params_zigzag(4:10:end,3),params_zigzag(5:10:end,3),params_zigzag(6:10:end,3),...
    params_zigzag(7:10:end,3),params_zigzag(8:10:end,3),params_zigzag(9:10:end,3),params_zigzag(10:10:end,3)];
values_ZZ.nugget = [params_zigzag(1:10:end,1),params_zigzag(1:10:end,4),params_zigzag(2:10:end,4), ...
    params_zigzag(3:10:end,4),params_zigzag(4:10:end,4),params_zigzag(5:10:end,4),params_zigzag(6:10:end,4),...
    params_zigzag(7:10:end,4),params_zigzag(8:10:end,4),params_zigzag(9:10:end,4),params_zigzag(10:10:end,4)];
values_ZZ.sill = [params_zigzag(1:10:end,1),params_zigzag(1:10:end,5),params_zigzag(2:10:end,5), ...
    params_zigzag(3:10:end,5),params_zigzag(4:10:end,5),params_zigzag(5:10:end,5),params_zigzag(6:10:end,5),...
    params_zigzag(7:10:end,5),params_zigzag(8:10:end,5),params_zigzag(9:10:end,5),params_zigzag(10:10:end,5)];
values_ZZ.all = params_zigzag;

clear c d data dim filename fig i j str T* zig* count_* opts lag_G vario* sill range nugget f myfit d_* params*
close all

%% Box of all swe options
clf
for g = 1:3
        glacier = char(options.glacier(g)); 
    cats    = {'S1','F1','S2','F2','S3','F3','S4','F4'};
    swedata = [];   group = [];     
    for opt = 2:9
        swedata = [swedata; sweOPT(opt).(glacier)(:,1)];
        group   = [group; repmat(cats(opt-1),length(sweOPT(opt).(glacier)),1)];
    end
    
%     [p,t,stats] = anova1(swedata, group);
%     multcompare(stats);
    
subplot(3,1,g)
    boxplot(swedata, group);
    ylim([0 1.2])
    ylabel('SWE (m w.e.)')
    
end

%Glacier 4
    textA = 'A                      A         A          A                      A         A';
    textB = '            B';
    textC = '                        C         C          C          C                     C';
        annotation('textbox',[0.16 0.885 0.1 0.1], 'String', textA,'EdgeColor','none')
        annotation('textbox',[0.16 0.870 0.1 0.1], 'String', textB,'EdgeColor','none')
        annotation('textbox',[0.16 0.855 0.1 0.1], 'String', textC,'EdgeColor','none')
        annotation('textbox',[0.01 0.855 0.1 0.1], 'String', '(a)','EdgeColor','none')
%Glacier 2
    textA = 'A          A          A                     A                      A         ';
    textB = 'B          B          B                     B          B          B         B';
    textC = '                                    C                     C                     C';
    textD = '            D                     D                      D                     D';
        annotation('textbox',[0.16 0.525 0.1 0.1], 'String', textA,'EdgeColor','none')
        annotation('textbox',[0.16 0.51 0.1 0.1], 'String', textB,'EdgeColor','none')
        annotation('textbox',[0.16 0.495 0.1 0.1], 'String', textC,'EdgeColor','none')
        annotation('textbox',[0.16 0.48 0.1 0.1], 'String', textD,'EdgeColor','none')
        annotation('textbox',[0.01 0.53 0.1 0.1], 'String', '(b)','EdgeColor','none')
%Glacier 13
    textA = 'A                     A                                            A         ';
    textB = '            B                      B                     B                    B';
    textC = 'C                     C                      C                     C        ';
    textD = '                       D                      D                     D';
        annotation('textbox',[0.16 0.225 0.1 0.1], 'String', textA,'EdgeColor','none')
        annotation('textbox',[0.16 0.21 0.1 0.1], 'String', textB,'EdgeColor','none')
        annotation('textbox',[0.16 0.195 0.1 0.1], 'String', textC,'EdgeColor','none')
        annotation('textbox',[0.16 0.18 0.1 0.1], 'String', textD,'EdgeColor','none')
        annotation('textbox',[0.01 0.23 0.1 0.1], 'String', '(c)','EdgeColor','none')
          

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 7 11];
saveFIG('AllSWEopts_boxplot')

clear text* fig filename cats swedata group g opt stats p t glacier


%% Variability for one measurement location

clf(1);
run MeasurementLocations.m  %This program determines the easting and northing of transect measurements
run Import_Density.m        %Imports snow density values
run Import_Transect.m       %Imports transect snow depth and measurement location data
glacier_list = ['G04';'G02';'G13']; %for selecting data from chosen glacier
for i = 1:3 %go through each glacier
    glacier = glacier_list(i,:);    
    G = char(options.glacier(i));
    z = pulldata(SD,'all',glacier,'all','all',1,'fat'); % pull transect data  

   data = z(5).depth(:,1:4);
   data = data*mean(cell2mat(Density.snowpit(:,2)))/1000; %SWE conversion
   data = (data - repmat(nanmean(data,2),1,4))./repmat(nanmean(data,2),1,4)*100;%./repmat(nanstd(data, [], 2),1,4);
         display([glacier,' 2 sigma ', num2str(2*nanstd(data(:)))]);
        %chi2gof(data(:))
   VARoneloc.(G) = data;
   
   
   if i==1;     bins = round(sqrt(length(SWEzz(i).swe)));   end
   edges   = linspace(-50,50,bins);
   N       =  histcounts(data,edges);
 figure(1);   plot((edges(:,1:end-1)+edges(:,2:end))/2,N/sum(N),'LineWidth',2, 'Color',options.RGB(i,:)); hold on 

stdtemp(i).swe = nanstd(z(5).depth(:,1:4), [], 2)./nanmean(z(5).depth(:,1:4),2)*100;
stdtemp(i).utm = z(5).depth(:,6:7);
end
            xlabel('SWE variability (%)'); ylabel('Frequency')
            legend('Glacier 4','Glacier 2','Glacier 13')
            %title('SWE (S1) variation at single measurement location')
            fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
            fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];
  saveFIG('SWEvarOneLocHIST')

 % Plot -> map of measurement std values
figure(2)
    topoParam.G4  = NaN(size(topo_full_ns.G4.elevation));
    topoParam.G2  = NaN(size(topo_full_ns.G2.elevation));
    topoParam.G13 = NaN(size(topo_full_ns.G13.elevation));

PlotTopoParameter(topoParam,'std in one grid cell', 'Coefficient of Variation (%)', stdtemp, 'colour', 'nomassB')
    saveFIG('Map_pointstd')

 
%% Varibaility in one DEM cell (multiple measurement locations)

clf
 %Boxplot
    T = [SWE(1).cellstd; SWE(2).cellstd; SWE(3).cellstd];
    G = [repmat('G04',length(SWE(1).cellstd),1); ...
         repmat('G02',length(SWE(2).cellstd),1);...
         repmat('G13',length(SWE(3).cellstd),1)];
boxplot(T,G,'Labels',{'Glacier 4','Glacier 2','Glacier 13'})
    ylabel([{'Standard deviation of SWE'},{'within DEM cell'}])
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',16)
     saveFIG('DEMcellSTD')
     
% PDF 
clf(1)
for g = 1:3;
glacier = char(options.glacier(g));

data = SWE(g).standard(~isnan(SWE(g).standard))./SWE(g).swe(~isnan(SWE(g).standard))*100;
        %chi2gof(data)
         display([glacier,' 2 sigma ', num2str(2*nanstd(data(:)))]);

if g == 1; bins = round(sqrt(numel(data))); end
   edges   = linspace(-30,30,bins);
   N       =  histcounts(data,edges);
 figure(1);   plot((edges(:,1:end-1)+edges(:,2:end))/2,N/sum(N),'LineWidth',2,'Color',options.RGB(g,:)); hold on 

stdtemp(g).swe = SWE(g).cellstd./SWE(g).swe*100;    stdtemp(g).utm = SWE(g).utm(:,1:2);

end
            xlabel('SWE variability (%)'); ylabel('Frequency')
            legend('Glacier 4','Glacier 2','Glacier 13')
            fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
            fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];    
            %title({'SWE variability due to multiple measurement','locations in one grid cell'})
saveFIG('SWEvarMeasureLocHIST')

% Plot -> map of cell std values
figure(2)
    topoParam.G4  = NaN(size(topo_full_ns.G4.elevation));
    topoParam.G2  = NaN(size(topo_full_ns.G2.elevation));
    topoParam.G13 = NaN(size(topo_full_ns.G13.elevation));

PlotTopoParameter(topoParam,'std in one grid cell', 'Coefficient of Variation (%)', stdtemp, 'colour', 'nomassB')
    saveFIG('Map_cellstd_measureLoc')
     
    
% Plot -> Hist of # of obs in one dem cell
figure(3)
for g = 1:3
   subplot(1,3,g)
   edges = 0.5:7.5;
   histogram(SWE(g).numObs,edges,'FaceColor',options.RGB(g,:));   
    xlim([0.5 7.5]); xlabel('Number of observations in a DEM cell')
    ylim([0 170]);                      ylabel('Count')
    title(options.glacier(g))
end
    saveFIG('NumObsPerCell')

%% Variability due to density options 

 %PDF
clf(1); 
for g = 1:3;
glacier = char(options.glacier(g));
clear data M S
    for i = 2:9
        data(:,i-1) = sweOPT(i).(glacier)(:,1);  
    end
 M = mean(data,2);      M = repmat(M,1,8);
 S = std(data,[],2);    S = repmat(S,1,8);
 data = (data-M)./M*100;%./S;
         display([glacier,' 2 sigma ', num2str(2*nanstd(data(:)))]);
 
if g ==1;  bins = round(sqrt(numel(data))/2); end
   edges   = linspace(-50,50,bins);
   N       =  histcounts(data,edges);
figure(1);   plot((edges(:,1:end-1)+edges(:,2:end))/2,N/sum(N),'LineWidth',2, 'Color',options.RGB(g,:)); hold on 
 
 % bins = round(sqrt(length(S)));
%     [N, edges] = histcounts(S(:,1),bins);
% figure(2);    plot((edges(:,1:end-1)+edges(:,2:end))/2,N,'LineWidth',2); hold on
%               xlabel('Standard Deviation due to Density Option'); ylabel('Frequency')
%               legend(options.glacier)
% mean(S(:,1))
end

figure(1); 
    xlabel('SWE variability (%)'); ylabel('Frequency')
    legend('Glacier 4','Glacier 2','Glacier 13') %title('SWE variation due to density interpolation')
            fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
            fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 6];    
saveFIG('SWEvarDensityHIST')

% for j = 1:length(data)
%     DD = data(:);%(data(j,:)-mean(data(j,:)))/std(data(j,:));
%     [N, edges] = histcounts(DD,bins);
%     plot((edges(:,1:end-1)+edges(:,2:end))/2,N,'LineWidth',2); hold on
% end

   
 %Map of std due to Density Interp
 figure(2)
for g = 1:3
    glacier = char(options.glacier(g));
    swedata.(glacier) = [];        
    for opt = 2:9
        swedata.(glacier) = [swedata.(glacier), sweOPT(opt).(glacier)(:,1)];
    end
    
stdtemp(g).swe = std(swedata.(glacier),[],2)./mean(swedata.(glacier),2)*100;    
stdtemp(g).utm = SWE(g).utm(:,1:2);    
end
    topoParam.G4  = NaN(size(topo_full_ns.G4.elevation));
    topoParam.G2  = NaN(size(topo_full_ns.G2.elevation));
    topoParam.G13 = NaN(size(topo_full_ns.G13.elevation));

PlotTopoParameter(topoParam,'std in one grid cell', 'Coefficient of Variation (%)', stdtemp, 'colour', 'nomassB')
    saveFIG('Map_cellstd_density')    