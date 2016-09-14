% Plotting Snow Density Data
%       This script used imported snow density data and plots density from
%       snowpit and SWE tube. Also plots against elevation

%       Inputs:         None
%       Other scripts:  SnowDensity_Import.m
%       Outputs:        None (just saves the plots)

%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Snowpit vs SWE tube density (all glaciers)
% errorbar(cell2mat(Density.pitANDtube(1:2,5)),cell2mat(Density.pitANDtube(1:2,2)), cell2mat(Density.pitANDtube(1:2,3)),'o'); hold on
% errorbar(cell2mat(Density.pitANDtube(3:4,5)),cell2mat(Density.pitANDtube(3:4,2)), cell2mat(Density.pitANDtube(3:4,3)),'o'); hold on
% errorbar(cell2mat(Density.pitANDtube(5:6,5)),cell2mat(Density.pitANDtube(5:6,2)), cell2mat(Density.pitANDtube(5:6,3)),'o'); hold on

    x = cell2mat(Density.pitANDtube(:,7)); %Snowpit
    y = cell2mat(Density.pitANDtube(:,2)); %SWE tube
    errory = [cell2mat(Density.pitANDtube(:,2))-cell2mat(Density.pitANDtube(:,4)), ...
        cell2mat(Density.pitANDtube(:,5))-cell2mat(Density.pitANDtube(:,2))]; %min and max tube
    errorx = [cell2mat(Density.pitANDtube(:,7))-cell2mat(Density.pitANDtube(:,9)), ...
        cell2mat(Density.pitANDtube(:,10))-cell2mat(Density.pitANDtube(:,7))]; %min and max SP
errorbarxy(x,y,errorx(:,1),errorx(:,2),errory(:,1),errory(:,2),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',5); hold on
    
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y);
plot(x,yfit,'r')
    xlabel('Snowpit density (kg m^{-3})')
    ylabel('SWE tube density (kg m^{-3})')
    dim = [0.15,0.5,0.11,0.11];
    str = {strcat('y= ',num2str(round(P(1),2)),'x + ',num2str(round(P(2),2))), ...
        strcat('R^2= ',num2str(round(LM.Rsquared.Ordinary,2)))} ;
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
    text(x+2, y+3, Density.pitANDtube(:,1))
    axis([290 400 220 400])
    axis equal
    
%filename = strcat('/home/glaciology1/Documents/Data/Plots/SnowpitVsSWEtube_all');
filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/SnowpitVsSWEtube_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error
    
%% SWEtube density vs elevation (all glaciers)
figure(2)
    x = nanmean(cell2mat(Density.tube(1:7,2:16)),2);
    y = cell2mat(Density.tube(1:7,22));
    errorx = nanstd(cell2mat(Density.tube(1:7,2:16)),1,2);
    errory = zeros(length(y),1);
h1 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','b','LineWidth',1,'MarkerSize',11); hold on

    x = nanmean(cell2mat(Density.tube(8:14,2:16)),2);
    y = cell2mat(Density.tube(8:14,22));
    errorx = nanstd(cell2mat(Density.tube(8:14,2:16)),1,2);
    errory = zeros(length(y),1);
h2 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',11); hold on

    x = nanmean(cell2mat(Density.tube(15:end,2:16)),2);
    y = cell2mat(Density.tube(15:end,22));
    errorx = nanstd(cell2mat(Density.tube(15:end,2:16)),1,2);
    errory = zeros(length(y),1);
h3 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','r','LineWidth',1,'MarkerSize',11); hold on

    xlabel('SWE tube density (kg m^{-3})')
    ylabel('Elevation(m)')
    legend([h1(1) h2(1) h3(1)], {'G04','G02','G13'})

%filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSWEtube_all');
filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/ElevationVsSWEtube_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error

%% Snowpit density vs elevation (all glaciers)
figure(3)
    x = cell2mat(Density.snowpit(1:4,2));
    y = cell2mat(Density.snowpit(1:4,5));
    errorx = 0.1*x;
    errory = zeros(length(y),1);
h2 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',11); hold on
    x = cell2mat(Density.snowpit(5:7,2));
    y = cell2mat(Density.snowpit(5:7,5));
    errorx = 0.1*x;
    errory = zeros(length(y),1);
h1 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','b','LineWidth',1,'MarkerSize',11); hold on   
    x = cell2mat(Density.snowpit(8:end,2));
    y = cell2mat(Density.snowpit(8:end,5));
    errorx = 0.1*x;
    errory = zeros(length(y),1);
h3 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','r','LineWidth',1,'MarkerSize',11); hold on

    xlabel('Snowpit density (kg m^{-3})')
    ylabel('Elevation(m)')
    legend([h1(1) h2(1) h3(1)], {'G04','G02','G13'})
%filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSnowpit_all');
filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/ElevationVsSnowpit_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error

%% Depth vs Density

% Probe depth
figure(1)
nonnanindex = [1:3,10:103];
[swefit, gof] = fit(cell2mat(Density.SWEdepth(nonnanindex,4)),cell2mat(Density.SWEdepth(nonnanindex,2)),'poly1');
plot(swefit, cell2mat(Density.SWEdepth(nonnanindex,4)),cell2mat(Density.SWEdepth(nonnanindex,2)));
    xlabel('SWE tube density (kg m^{-3})')
    ylabel('Mean snow probe depth (cm)')
    title({'Depth vs Density','(Snow probe depth)'});
    dim = [.20 .5 .3 .3];
    str = ['R^2 = ', num2str(round(gof.rsquare,2))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
    p = legend('data');
    set(p,'visible','off')
    
% SWE tube depth
figure(2)
nonnanindex = [1:103];
[swefit, gof] = fit(cell2mat(Density.SWEdepth(nonnanindex,4)),cell2mat(Density.SWEdepth(nonnanindex,3)),'poly1');
plot(swefit, cell2mat(Density.SWEdepth(nonnanindex,4)),cell2mat(Density.SWEdepth(nonnanindex,3)));
    xlabel('SWE tube density (kg m^{-3})')
    ylabel('SWE tube depth (cm)')
    title({'Depth vs Density','(SWE tube depth)'});
    dim = [.20 .5 .3 .3];
    str = ['R^2 = ', num2str(round(gof.rsquare,2))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
    p = legend('data');
    set(p,'visible','off')    
    
% Snowpit
figure(3)
[swefit, gof] = fit(cell2mat(Density.snowpit(:,2)),cell2mat(Density.snowpit(:,8)),'poly1');
plot(swefit, cell2mat(Density.snowpit(:,2)),cell2mat(Density.snowpit(:,8)));
    xlabel('Snowpit density (kg m^{-3})')
    ylabel('Snowpit depth (cm)')
    title({'Depth vs Density','(Snowpit)'});
    dim = [.20 .5 .3 .3];
    str = ['R^2 = ', num2str(round(gof.rsquare,2))];
    annotation('textbox',dim,'String',str,'FitBoxToText','on')
    p = legend('data');
    set(p,'visible','off')      
    
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
display(['mean = ',num2str(round(nanmean(nanmean(cell2mat(Density.tube(range,2:19)))))),...
        ', std = ',num2str(round(nanstd(nanstd(cell2mat(Density.tube(range,2:19)))))),...
        ', n = ',num2str(length(range))]);
display('G02 tube'); range = 8:14;
display(['mean = ',num2str(round(nanmean(nanmean(cell2mat(Density.tube(range,2:19)))))),...
        ', std = ',num2str(round(nanstd(nanstd(cell2mat(Density.tube(range,2:19)))))),...
        ', n = ',num2str(length(range))]);
display('G13 tube'); range = 15:31;
display(['mean = ',num2str(round(nanmean(nanmean(cell2mat(Density.tube(range,2:19)))))),...
        ', std = ',num2str(round(nanstd(nanstd(cell2mat(Density.tube(range,2:19)))))),...
        ', n = ',num2str(length(range))]);
display('All tube'); range = 1:31;
display(['mean = ',num2str(round(nanmean(nanmean(cell2mat(Density.tube(range,2:19)))))),...
        ', std = ',num2str(round(nanstd(nanstd(cell2mat(Density.tube(range,2:19)))))),...
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
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2:19)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G02 tube'); range = 8:14;
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2:19)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G13 tube'); range = [15,17,18,20:33];
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2:19)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
%% SWE depth correction    

    depth = cell2mat(Density.SWEdepth(:,3));
    den = cell2mat(Density.SWEdepth(:,4));
C = [depth, den];
[sortC, I] = sortrows(C);
detrendC = [I, detrend(C(:,2)) + mean(C(:,2))];
final = sortrows(detrendC);
%plot(final(:,2),'.')

tube = [Density.SWEdepth(4,1), mean(final(4:9,2)), std(final(4:9,2)); Density.SWEdepth(13,1), mean(final(13:19,2)), std(final(13:19,2));...
    Density.SWEdepth(29,1), mean(final(29:35,2)), std(final(29:35,2)); Density.SWEdepth(44,1), mean(final(44:50,2)), std(final(44:50,2));...
    Density.SWEdepth(54,1), mean(final(54:61,2)), std(final(54:61,2)); Density.SWEdepth(85,1), mean(final(84:89,2)), std(final(84:89,2))];
tube(:,2:3) = num2cell(round(cell2mat(tube(:,2:3))));

tubeG = [cellstr('G04'), mean([final(4:9,2);final(13:19,2)]), std([final(4:9,2);final(13:19,2)]);...
    cellstr('G02'), mean([final(29:35,2);final(44:50,2)]), std([final(29:35,2);final(44:50,2)]);...
    cellstr('G13'), mean([final(54:61,2);final(84:89,2)]), std([final(54:61,2);final(84:89,2)]);...
    cellstr('All'), mean([final(4:9,2);final(13:19,2);final(29:35,2);final(44:50,2);final(54:61,2);final(84:89,2)]),...
        std([final(4:9,2);final(13:19,2);final(29:35,2);final(44:50,2);final(54:61,2);final(84:89,2)])];
tubeG(:,2:3) = num2cell(round(cell2mat(tubeG(:,2:3))));

    x = cell2mat(Density.pitANDtube(:,7)); %Snowpit
    y = cell2mat(tube(:,2)); %Detrended SWE tube
    errory = [0, 0];
    %errory = [cell2mat(Density.pitANDtube(:,2))-cell2mat(Density.pitANDtube(:,4)), ...
    %    cell2mat(Density.pitANDtube(:,5))-cell2mat(Density.pitANDtube(:,2))]; %min and max tube
    errorx = [cell2mat(Density.pitANDtube(:,7))-cell2mat(Density.pitANDtube(:,9)), ...
        cell2mat(Density.pitANDtube(:,10))-cell2mat(Density.pitANDtube(:,7))]; %min and max SP
errorbarxy(x,y,errorx(:,1),errorx(:,2),errory(:,1),errory(:,2),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',5); hold on
    
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y);
plot(x,yfit,'r')
    xlabel('Snowpit density (kg m^{-3})')
    ylabel('SWE tube density (kg m^{-3})')
    dim = [0.15,0.5,0.11,0.11];
    str = {strcat('y=',num2str(round(P(1),2)),'*x+',num2str(round(P(2),2))), ...
        strcat('R^2=',num2str(round(LM.Rsquared.Ordinary,2)))} ;
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
    text(x+2, y+2, Density.pitANDtube(:,1))
%     axis([290 400 220 400])
    axis equal

%% Density Range

%Snowpit
    %Name, ref, min ,max, range
disp = [Density.snowpit(:,1), num2cell(round(cell2mat(Density.snowpit(:,2)))),num2cell(round(cell2mat(Density.snowpit(:,6)))), num2cell(round(cell2mat(Density.snowpit(:,7))))];
disp = [disp, num2cell(cell2mat(disp(:,4))-cell2mat(disp(:,3)))];

%Tube
    %Name, num, mean density, min, max, range as % of mean
disp2 = [Density.tube(:,1), num2cell(sum(~isnan(cell2mat(Density.tube(:,2:18))),2)), num2cell(round(nanmean(cell2mat(Density.tube(:,2:18)),2))),...
    num2cell(round(nanmin(cell2mat(Density.tube(:,2:18)),[],2))),num2cell(round(nanmax(cell2mat(Density.tube(:,2:18)),[],2)))];
disp2(:,6) = num2cell(round((cell2mat(disp2(:,5))-cell2mat(disp2(:,4)))./cell2mat(disp2(:,3))*100));

