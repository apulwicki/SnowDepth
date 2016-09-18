%% Point Scale

% Box plots for each glacier

boxplot([SWE(1).depth;SWE(2).depth;SWE(3).depth], [SWE(1).glacier;SWE(2).glacier;SWE(3).glacier],...
    'GroupOrder',{'G04','G02','G13'},'labels',{'Glacier 4','Glacier 2','Glacier 13'})
    ylabel('Snow depth (cm)')
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',12) 
    
    print([options.path1,'box_depth'],'-dpng'); print([options.path2,'box_depth'],'-dpng') 

%% Zigzag

zig_lab = ['G04\_Z3A\_ZZ0'; 'G04\_Z2A\_ZZ0'; 'G04\_Z5B\_ZZ0';...
            'G02\_Z5C\_ZZ0'; 'G02\_Z7A\_ZZ0'; 'G02\_Z3B\_ZZ0'; ...
            'G13\_Z7C\_ZZ0'; 'G13\_Z4C\_ZZ0'; 'G13\_Z3B\_ZZ0'; 'G13\_Z5A\_ZZ0'];
for j = 1:3
    T1 = cellstr(char(SWE(j).label));

    for i = 1:length(zig_lab)
        T2 = find(~cellfun('isempty',strfind(T1,zig_lab(i,:))));
        if isempty(T2)
            break
        end  

        data = [SWE(j).swe(T2)*100, SWE(j).utm(T2,1:2)];
            data(:,2) = max(data(:,2))-data(:,2);
            data(:,3) = max(data(:,3))-data(:,3);

        scatter(data(:,2),data(:,3),25, data(:,1),'filled')
            title(zig_lab(i,1:8))
            str = {strcat('mean= ', num2str(round(mean(data(:,1)),2)),'cm SWE'), ...
                strcat('std= ', num2str(round(std(data(:,1)),2)),'cm SWE')};
            dim = [.13 .5 .3 .3];
            annotation('textbox',dim,'String', str,'FitBoxToText','on')
            xlabel('Distance (m)')
            ylabel('Distance (m)')
                c = colorbar;
                c.Label.String = 'SWE (cm)';
            filename = [zig_lab(i,1:3), zig_lab(i,6:8)];
            print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
            clf

        d = variogramAlex(data, 2, 40);
        fit = variofitAlex(d,zig_lab(i,1:8));

            filename = [filename, 'variogram'];
            fig = gcf; fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 8 9];
            print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
            clf

    end
end