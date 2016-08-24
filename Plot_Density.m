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

    x = cell2mat(Density.pitANDtube(:,5));
    y = cell2mat(Density.pitANDtube(:,2));
    errory = cell2mat(Density.pitANDtube(:,3));
    errorx = .1*x;
errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',5); hold on
    
    x = cell2mat(Density.pitANDtube(:,5)); y = cell2mat(Density.pitANDtube(:,2));
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y);
plot(x,yfit,'r')
    xlabel('Snowpit density (kg m^{-3})')
    ylabel('SWE tube density (kg m^{-3})')
    dim = [0.15,0.5,0.11,0.11];
    str = {strcat('y=',num2str(round(P(1),2)),'*x+',num2str(round(P(2),2))), ...
        strcat('R^2=',num2str(round(LM.Rsquared.Ordinary,2)))} ;
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
    text(cell2mat(Density.pitANDtube(:,5))+3, cell2mat(Density.pitANDtube(:,2))+2, Density.pitANDtube(:,1))
    axis([230 400 230 400])
    
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
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',11); hold on

    x = nanmean(cell2mat(Density.tube(8:14,2:16)),2);
    y = cell2mat(Density.tube(8:14,22));
    errorx = nanstd(cell2mat(Density.tube(8:14,2:16)),1,2);
    errory = zeros(length(y),1);
h2 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','b','LineWidth',1,'MarkerSize',11); hold on

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
h1 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',11); hold on
    x = cell2mat(Density.snowpit(5:7,2));
    y = cell2mat(Density.snowpit(5:7,5));
    errorx = 0.1*x;
    errory = zeros(length(y),1);
h2 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
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
