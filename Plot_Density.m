% Plotting Snow Density Data
%       This script used imported snow density data and plots density from
%       snowpit and SWE tube. Also plots against elevation

%       Inputs:         None
%       Other scripts:  SnowDensity_Import.m
%       Outputs:        None (just saves the plots)

%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Snowpit vs SWE tube density (all glaciers)
errorbar(cell2mat(Density.pitANDtube(1:2,5)),cell2mat(Density.pitANDtube(1:2,2)), cell2mat(Density.pitANDtube(1:2,3)),'o'); hold on
errorbar(cell2mat(Density.pitANDtube(3:4,5)),cell2mat(Density.pitANDtube(3:4,2)), cell2mat(Density.pitANDtube(3:4,3)),'o'); hold on
errorbar(cell2mat(Density.pitANDtube(5:6,5)),cell2mat(Density.pitANDtube(5:6,2)), cell2mat(Density.pitANDtube(5:6,3)),'o'); hold on
    x = cell2mat(Density.pitANDtube(:,5)); y = cell2mat(Density.pitANDtube(:,2));
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y);
plot(x,yfit,'r')
    xlabel('Snowpit density (kg m^{-3})')
    ylabel('SWE tube density (kg m^{-3})')
    legend('G04','G02','G13','Fit')
    dim = [0.15,0.5,0.11,0.11];
    str = {strcat('y=',num2str(round(P(1),2)),'*x+',num2str(round(P(2),2))), ...
        strcat('R^2=',num2str(round(LM.Rsquared.Ordinary,2)))} ;
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
    text(cell2mat(Density.pitANDtube(:,5))+1, cell2mat(Density.pitANDtube(:,2))-1, Density.pitANDtube(:,1))
    
filename = strcat('/home/glaciology1/Documents/Data/Plots/SnowpitVsSWEtube_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error

%% SWEtube density vs elevation (all glaciers)
figure(2)
    x = nanmean(cell2mat(Density.tube(1:7,2:16)),2);
    y = cell2mat(Density.tube(1:7,22));
    error = nanstd(cell2mat(Density.tube(1:7,2:16)),1,2);
errorbarxy(x,y,error,[],{'ko', 'b', 'r'}); hold on

    x = nanmean(cell2mat(Density.tube(8:14,2:16)),2);
    y = cell2mat(Density.tube(8:14,22));
    error = nanstd(cell2mat(Density.tube(8:14,2:16)),1,2);
errorbarxy(x,y,error,[],{'ro', 'b', 'r'}); hold on

    x = nanmean(cell2mat(Density.tube(15:end,2:16)),2);
    y = cell2mat(Density.tube(15:end,22));
    error = nanstd(cell2mat(Density.tube(15:end,2:16)),1,2);
errorbarxy(x,y,error,[],{'mo', 'b', 'r'}); hold on

    xlabel('SWE tube density (kg m^{-3})')
    ylabel('Elevation(m)')
    legend('G04','G02','G13')
filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSWEtube_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error

%% Snowpit density vs elevation (all glaciers)
figure(3)
    x = cell2mat(Density.snowpit(1:4,2));
    y = cell2mat(Density.snowpit(1:4,5));
scatter(x,y,'Filled'); hold on
    x = cell2mat(Density.snowpit(5:7,2));
    y = cell2mat(Density.snowpit(5:7,5));
scatter(x,y,'Filled'); hold on    
    x = cell2mat(Density.snowpit(8:end,2));
    y = cell2mat(Density.snowpit(8:end,5));
scatter(x,y,'Filled')

    xlabel('Snowpit density (kg m^{-3})')
    ylabel('Elevation(m)')
    legend('G02','G04','G13')
filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSnowpit_all');
print(filename,'-dpng')
    clear P LM x y index* i j count yfit str dim filename error
