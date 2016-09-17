%% Point Scale

% Box plots for each glacier

boxplot([SWE(1).depth;SWE(2).depth;SWE(3).depth], [SWE(1).glacier;SWE(2).glacier;SWE(3).glacier],...
    'GroupOrder',{'G04','G02','G13'},'labels',{'Glacier 4','Glacier 2','Glacier 13'})
    ylabel('Snow depth (cm)')
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',12) 
    
    filename = strcat('/home/glaciology1/Documents/Data/Plots/box_depth');
    %filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/box_depth');
    print(filename,'-dpng')
    filename = strcat('/home/glaciology1/Documents/MastersDocuments/Methods/box_depth');
    %filename = strcat('/Users/Alexandra/Documents/SFU/MastersDocuments/Methods/box_depth');
    print(filename,'-dpng')

%% Zigzag

    %pulldataSWE(data, glacier, pattern, book, person)
z = pulldataSWE(SWE, 'all','all','all','all');
    %scatter(z.utm(:,1),z.utm(:,2), 10, z.depth,'filled');

d = variogramAlex([z x1 y1], 15, 'default');
fit = variofitAlex(d,glacier);

    filename = strcat('/home/glaciology1/Documents/Data/Plots/variofull',glacier);
   % filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/variofull',glacier);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 9];
    print(filename,'-dpng','-r0')
    clf