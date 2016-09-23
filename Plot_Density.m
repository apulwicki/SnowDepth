% Plotting Snow Density Data
%       This script used imported snow density data and plots density from
%       snowpit and SWE tube. Also plots against elevation

%       Inputs:         None
%       Other scripts:  SnowDensity_Import.m
%       Outputs:        None (just saves the plots)

%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Snowpit vs SWE tube density (all glaciers)

    x = cell2mat(Density.pitANDtube(:,7)); %Snowpit
    y = cell2mat(Density.pitANDtube(:,2)); %SWE tube
    errory = [cell2mat(Density.pitANDtube(:,2))-cell2mat(Density.pitANDtube(:,4)), ...
        cell2mat(Density.pitANDtube(:,5))-cell2mat(Density.pitANDtube(:,2))]; %min and max tube
    errorx = [cell2mat(Density.pitANDtube(:,7))-cell2mat(Density.pitANDtube(:,9)), ...
        cell2mat(Density.pitANDtube(:,10))-cell2mat(Density.pitANDtube(:,7))]; %min and max SP
errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',5); hold on
    
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y);
plot(x,yfit,'r')
    xlabel('Snowpit density (kg m^{-3})')
    ylabel('SWE tube density (kg m^{-3})')
    dim = [0.15,0.16,0.11,0.11];
    str = {strcat('y= ',num2str(round(P(1),2)),'x+ ',num2str(round(P(2)))), ...
        strcat('R^2= ',num2str(round(LM.Rsquared.Ordinary,2)))} ;
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
    text(x+2, y+3, Density.pitANDtube(:,1))
    axis([290 400 220 400])
    axis equal
    
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',12) 

    filename = 'SnowpitVsSWEtube_all';
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')


    clear P LM x y index* i j count yfit str dim filename error* ans
    
%% SWEtube density vs elevation (all glaciers)

    y = cell2mat(Density.tube(1:7,2)); %Glacier 4
    x = cell2mat(Density.tube(1:7,9));
    errory = [cell2mat(Density.tube(1:7,2))-cell2mat(Density.tube(1:7,4)),...
        cell2mat(Density.tube(1:7,5))-cell2mat(Density.tube(1:7,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit1, gof1] = fit(x,y,'poly1');
h1 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','b','LineWidth',1,'MarkerSize',11); hold on
l1 = plot(SPfit1,'b'); hold on

    y = cell2mat(Density.tube(8:14,2)); %Glacier 2
    x = cell2mat(Density.tube(8:14,9));
    errory = [cell2mat(Density.tube(8:14,2))-cell2mat(Density.tube(8:14,4)),...
        cell2mat(Density.tube(8:14,5))-cell2mat(Density.tube(8:14,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit2, gof2] = fit(x,y,'poly1');
h2 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',11); hold on
l2 = plot(SPfit2,'k'); hold on

    y = cell2mat(Density.tube(15:end,2)); %Glacier 13
    x = cell2mat(Density.tube(15:end,9));
    errory = [cell2mat(Density.tube(15:end,2))-cell2mat(Density.tube(15:end,4)),...
        cell2mat(Density.tube(15:end,5))-cell2mat(Density.tube(15:end,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit3, gof3] = fit(x,y,'poly1');
h3 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','r','LineWidth',1,'MarkerSize',11); hold on
l3 = plot(SPfit3,'r'); hold on

    ylabel('SWE tube density (kg m^{-3})')
    xlabel('Elevation(m)')
    legend([h1(1) l1 h2(1) l2 h3(1) l3], {'G04',['G04 fit R^2=',num2str(round(gof1.rsquare,2))], ...
        'G02',['G02 fit R^2=',num2str(round(gof2.rsquare,2))],...
        'G13',['G13 fit R^2=',num2str(round(gof3.rsquare,2))]})
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',12) 
    
filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSWEtube_all');
%filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/ElevationVsSWEtube_all');
print(filename,'-dpng')
filename = strcat('/home/glaciology1/Documents/MastersDocuments/Methods/ElevationVsSWEtube_all');
%filename = strcat('/Users/Alexandra/Documents/SFU/MastersDocuments/Methods/ElevationVsSWEtube_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error

    
display('G04 SP'); 
display([num2str(round(SPfit1.p1,3)),'x + ',num2str(round(SPfit1.p2,1))])

display('G02 SP'); 
display([num2str(round(SPfit2.p1,3)),'x + ',num2str(round(SPfit2.p2,1))])

display('G13 SP');
display([num2str(round(SPfit3.p1,3)),'x + ',num2str(round(SPfit3.p2,1))])
    
%% Snowpit density vs elevation (all glaciers)
figure(3)
    y = cell2mat(Density.snowpit(1:4,2));
    x = cell2mat(Density.snowpit(1:4,5));
    errory = [cell2mat(Density.snowpit(1:4,2))-cell2mat(Density.snowpit(1:4,6)),...
        cell2mat(Density.snowpit(1:4,7))-cell2mat(Density.snowpit(1:4,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit2, gof2] = fit(x,y,'poly1');
h2 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',11); hold on
l2 = plot(SPfit2,'k'); hold on

    y = cell2mat(Density.snowpit(5:7,2));
    x = cell2mat(Density.snowpit(5:7,5));
    errory = [cell2mat(Density.snowpit(5:7,2))-cell2mat(Density.snowpit(5:7,6)),...
        cell2mat(Density.snowpit(5:7,7))-cell2mat(Density.snowpit(5:7,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit1, gof1] = fit(x,y,'poly1');
h1 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','b','LineWidth',1,'MarkerSize',11); hold on   
l1 = plot(SPfit1,'b'); hold on

    y = cell2mat(Density.snowpit(8:end,2));
    x = cell2mat(Density.snowpit(8:end,5));
    errory = [cell2mat(Density.snowpit(8:end,2))-cell2mat(Density.snowpit(8:end,6)),...
        cell2mat(Density.snowpit(8:end,7))-cell2mat(Density.snowpit(8:end,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit3, gof3] = fit(x,y,'poly1');
h3 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','r','LineWidth',1,'MarkerSize',11); hold on
l3 = plot(SPfit3,'r'); hold on

    ylabel('Snowpit density (kg m^{-3})')
    xlabel('Elevation(m)')
    legend([h1(1) l1 h2(1) l2 h3(1) l3], {'G04',['G04 fit R^2=',num2str(round(gof1.rsquare,2))], ...
            'G02',['G02 fit R^2=',num2str(round(gof2.rsquare,2))],...
            'G13',['G13 fit R^2=',num2str(round(gof3.rsquare,2))]},'Location','best')
filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSnowpit_all');
%filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/ElevationVsSnowpit_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error

    
display('G04 SP'); 
display([num2str(round(SPfit1.p1,3)),'x + ',num2str(round(SPfit1.p2,1))])
    
display('G02 SP'); 
display([num2str(round(SPfit2.p1,3)),'x + ',num2str(round(SPfit2.p2,1))])

display('G13 SP');
display([num2str(round(SPfit3.p1,3)),'x + ',num2str(round(SPfit3.p2,1))])
%% Depth vs Density

figure(2)
nonnanindex = [1:103];
[swefit, gof] = fit(cell2mat(Density.SWEdepth(nonnanindex,4)),cell2mat(Density.SWEdepth(nonnanindex,3)),'poly1');
plot(swefit); hold on
plot(cell2mat(Density.SWEdepth(nonnanindex,4)),cell2mat(Density.SWEdepth(nonnanindex,3)),'.','MarkerSize',13);
    xlabel('SWE tube density (kg m^{-3})')
    ylabel('SWE tube depth (cm)')
    %title({'Depth vs Density','(SWE tube depth)'});
    dim = [.20 .5 .3 .3];
    str = ['R^2 = ', num2str(round(gof.rsquare,2))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
    p = legend('data');
    set(p,'visible','off')    
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',12) 
filename = strcat('/home/glaciology1/Documents/Data/Plots/DepthDensity_SWEtube');
%filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/DepthDensity_SWEtube');
print(filename,'-dpng')
filename = strcat('/home/glaciology1/Documents/MastersDocuments/Methods/DepthDensity_SWEtube');
%filename = strcat('/Users/Alexandra/Documents/SFU/MastersDocuments/Methods/DepthDensity_SWEtube');
print(filename,'-dpng')

% Snowpit
figure(3)
[swefit, gof] = fit(cell2mat(Density.snowpit(:,2)),cell2mat(Density.snowpit(:,8)),'poly1');
plot(swefit); hold on
plot(cell2mat(Density.snowpit(:,2)),cell2mat(Density.snowpit(:,8)),'.','markers',13);
    xlabel('Snowpit density (kg m^{-3})')
    ylabel('Snowpit depth (cm)')
    %title({'Depth vs Density','(Snowpit)'});
    dim = [.20 .5 .3 .3];
    str = ['R^2 = ', num2str(round(gof.rsquare,2))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
    p = legend('data');
    set(p,'visible','off')      
    fig=gcf;
    set(findall(fig,'-property','FontSize'),'FontSize',12) 
filename = strcat('/home/glaciology1/Documents/Data/Plots/DepthDensity_SP');
%filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/DepthDensity_SP');
print(filename,'-dpng')      
filename = strcat('/home/glaciology1/Documents/MastersDocuments/Methods/DepthDensity_SP');
%filename = strcat('/Users/Alexandra/Documents/SFU/MastersDocuments/Methods/DepthDensity_SP');
print(filename,'-dpng')


%% Depth vs Elevation
% Swe tube
depth = cell2mat(Density.tube(:,10));
elev = cell2mat(Density.tube(:,9));
index = [1,7; 8,14; 15,31];
for i = 1:3
   plot(depth(index(i,1):index(i,2)), elev(index(i,1):index(i,2)), '.','MarkerSize', 15); hold on
end
    xlabel('Depth (cm)'); ylabel('Elevation (m a.s.l.)');
    legend('Glacier 4','Glacier 2','Glacier 13')
    
% All data
colour = ['b','k','r'];
for j = 1:3
    depth = SWE(j).depth(~isnan(SWE(j).utm(:,3)));
    elev = SWE(j).utm(~isnan(SWE(j).utm(:,3)),3);
    i = ['n',num2str(j)];
    [f.(i), g.(i)] = fit(depth,elev,'poly1');
        mark = ['.',colour(j)];
    h = plot(depth,elev,mark,'MarkerSize',13); hold on
    plot(f.(i),colour(j)); hold on
end
    xlabel('Depth (cm)'); ylabel('Elevation (m a.s.l.)');
    legend('Glacier 4',['R^2=',num2str(round(g.n1.rsquare,2))],...
        'Glacier 2',['R^2=',num2str(round(g.n2.rsquare,2))],...
        'Glacier 13',['R^2=',num2str(round(g.n3.rsquare,2))],'Location','best')
%% Basic stats

% Snowpit
clc
display('G04 SP'); range = 5:7;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G02 SP'); range = 1:4;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G13 SP'); range = 8:10;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);
display('All SP'); range = 1:10;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);

% SWE tube
display('')
display('G04 tube'); range = 1:7;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G02 tube'); range = 8:14;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G13 tube'); range = 15:31;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
display('All tube'); range = 1:31;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
    
%% Linear regression
clc

%SP
display('G04 SP'); range = 5:7;
    SPfit = fit(cell2mat(Density.snowpit(range,5)), cell2mat(Density.snowpit(range,2)),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G02 SP'); range = 1:4;
    SPfit = fit(cell2mat(Density.snowpit(range,5)), cell2mat(Density.snowpit(range,2)),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])

display('G13 SP'); range = 8:10;
    SPfit = fit(cell2mat(Density.snowpit(range,5)), cell2mat(Density.snowpit(range,2)),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])

%SP
display('')
display('G04 tube'); range = 1:7;
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G02 tube'); range = 8:14;
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G13 tube'); range = [15,17,18,20:33];
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])

%% Density Range

%Snowpit
    %Name, ref, min ,max, range
disp = [Density.snowpit(:,1), num2cell(round(cell2mat(Density.snowpit(:,2)))),num2cell(round(cell2mat(Density.snowpit(:,6)))), num2cell(round(cell2mat(Density.snowpit(:,7))))];
disp = [disp, num2cell(cell2mat(disp(:,4))-cell2mat(disp(:,3)))];

%Tube
disp2 = [Density.tube(:,1), Density.tube(:,6), Density.tube(:,2), Density.tube(:,4:5)...
    num2cell(round((cell2mat(Density.tube(:,5))-cell2mat(Density.tube(:,4)))/cell2mat(Density.tube(:,2))*100))];
disp2(:,6:11) = [];