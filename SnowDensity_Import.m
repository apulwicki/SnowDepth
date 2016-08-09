% Snow density import

[snowpit_density, snowpit_text, snowpit_raw] = xlsread('summary_densitydata.xlsx','snowpit','A2:E11');

[SWEtube_density, SWEtube_text, SWEtube_raw] = xlsread('summary_densitydata.xlsx','SWEtube','A2:V34');

%Establishing index for various data
index_zigzag = [1,3,5,8,10,13,26,28,30,32];
index_SPSWE = [2,4,9,14,15,27];

%Select only good qulaity data
for i = 1:size(SWEtube_density,1)
    for j = 1:size(SWEtube_density,2)
        if SWEtube_density(i,j)==0
           SWEtube_density(i,j-1:j)=nan; 
        elseif SWEtube_density(i,j)==1
           SWEtube_density(i,j)=nan; 
        end
    end
end


%% Zigzag SWE locations

zigzagSWE = SWEtube_density(index_zigzag',:);

    %Select only good qulaity data
count = zeros(size(zigzagSWE,1),1);
for i = 1:size(zigzagSWE,1)
    count(i,1) = sum(~isnan(zigzagSWE(i,1:15)));
end

    %Get mean density and std
    % Location, mean density, std, number of obs
zigzagSWE = [SWEtube_text(index_zigzag',:), num2cell(nanmean(zigzagSWE(:,1:15),2)/1000),...
    num2cell(nanstd(zigzagSWE(:,1:15),1,2)/1000), num2cell(count)];

clear i j count

%% Snowpit vs SWE tube

%Snowpit SWE tube values
snowpittubeSWE = SWEtube_density(index_SPSWE',:);

    %Select only good qulaity data
count = zeros(size(snowpittubeSWE,1),1);
for i = 1:size(snowpittubeSWE,1)
    for j = 1:size(snowpittubeSWE,2)
        if snowpittubeSWE(i,j)==0
           snowpittubeSWE(i,j-1:j)=nan; 
        elseif snowpittubeSWE(i,j)==1
           snowpittubeSWE(i,j)=nan; 
        end
    end
    count(i,1) = sum(~isnan(snowpittubeSWE(i,1:15)));
end

    %Get mean density and std
    % Location, mean density, std, number of obs, snowpit density
snowpittubeSWE = [SWEtube_text(index_SPSWE',:), num2cell(nanmean(snowpittubeSWE(:,1:15),2)/1000),...
    num2cell(nanstd(snowpittubeSWE(:,1:15),1,2)/1000), num2cell(count)];

    %Corresponding snowpit densities
indexSP = zeros(size(index_SPSWE,2),1);
for i = 1:size(index_SPSWE,2)
    indexSP(i,1) = find(strcmp(snowpittubeSWE(i,1),snowpit_text));
end
snowpittubeSWE = [snowpittubeSWE, num2cell(snowpit_density(indexSP,1)/1000)];

%% Plots

%Snowpit vs SWE tube (all G)
errorbar(cell2mat(snowpittubeSWE(1:2,5)),cell2mat(snowpittubeSWE(1:2,2)), cell2mat(snowpittubeSWE(1:2,3)),'o'); hold on
errorbar(cell2mat(snowpittubeSWE(3:4,5)),cell2mat(snowpittubeSWE(3:4,2)), cell2mat(snowpittubeSWE(3:4,3)),'o'); hold on
errorbar(cell2mat(snowpittubeSWE(5:6,5)),cell2mat(snowpittubeSWE(5:6,2)), cell2mat(snowpittubeSWE(5:6,3)),'o'); hold on
    x = cell2mat(snowpittubeSWE(:,5)); y = cell2mat(snowpittubeSWE(:,2));
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y);
plot(x,yfit,'r')
    xlabel('Snowpit density (kg^3 m^{-3})')
    ylabel('SWE tube density (kg^3 m^{-3})')
    legend('G04','G02','G13','Fit')
    dim = [0.15,0.5,0.11,0.11];
    str = {strcat('y=',num2str(round(P(1),2)),'*x+',num2str(round(P(2),2))), ...
        strcat('R^2=',num2str(round(LM.Rsquared.Ordinary,2)))} ;
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
    
filename = strcat('/home/glaciology1/Documents/Data/Plots/SnowpitVsSWEtube_all');
print(filename,'-dpng')

%SWEtube vs elevation (all G)
figure(2)
    x = nanmean(SWEtube_density(1:7,1:15),2);
    y = SWEtube_density(1:7,21);
    error = nanstd(SWEtube_density(1:7,1:15),1,2);
errorbarxy(x,y,error,[],{'ko', 'b', 'r'}); hold on

    x = nanmean(SWEtube_density(8:14,1:15),2);
    y = SWEtube_density(8:14,21);
    error = nanstd(SWEtube_density(8:14,1:15),1,2);
errorbarxy(x,y,error,[],{'ro', 'b', 'r'}); hold on

    x = nanmean(SWEtube_density(15:end,1:15),2);
    y = SWEtube_density(15:end,21);
    error = nanstd(SWEtube_density(15:end,1:15),1,2);
errorbarxy(x,y,error,[],{'mo', 'b', 'r'}); hold on

    xlabel('SWE tube density (kg m^{-3})')
    ylabel('Elevation(m)')
    legend('G04','G02','G13')
filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSWEtube_all');
print(filename,'-dpng')
%%
%Snowpit vs elevation (all G)
figure(3)
    x = snowpit_density(1:4,1);
    y = snowpit_density(1:4,4);
scatter(x,y,'Filled'); hold on
    x = snowpit_density(5:7,1);
    y = snowpit_density(5:7,4);
scatter(x,y,'Filled'); hold on    
    x = snowpit_density(8:end,1);
    y = snowpit_density(8:end,4);
scatter(x,y,'Filled')

    xlabel('Snowpit density (kg m^{-3})')
    ylabel('Elevation(m)')
    legend('G02','G04','G13')
filename = strcat('/home/glaciology1/Documents/Data/Plots/ElevationVsSnowpit_all');
print(filename,'-dpng')


clear P LM x y index* i j count yfit str dim filename
%%
% Clear variables that you don't want exported
%clear SWEtube_* snowpit_*
