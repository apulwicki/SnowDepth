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


%% Std in one DEM cell

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