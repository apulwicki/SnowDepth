% Importing Snow Density Data
%       This script imports snow density data from the summary of density
%       data file. It also differentes the SWE tube density values and the
%       snowpit density values.

%       Inputs:         Summary of density data ('summary_densitydata.xlsx')
%       Outputs:        Density data structure with separated values for
%                           just snowpit, snowpit and tube, just tube, and
%                           tube in zigzag locations. 

%       Alexandra Pulwicki  Created: August 2016
%                           Updated: September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing density values

%Get data from file
    [snowpit_density, snowpit_text, snowpit_raw] = xlsread('summary_densitydata.xlsx','snowpit','A2:H11');
    [SWEtube_density, SWEtube_text, SWEtube_raw] = xlsread('summary_densitydata.xlsx','SWEtube','A2:V34');
    [depth_density, depth_text, depth_raw] = xlsread('summary_densitydata.xlsx','DepthDensity','A2:D107');
    
%Select only good quality data
    for i = 1:size(SWEtube_density,1) %for all elemnts in the matrix
        for j = 1:size(SWEtube_density,2)
            if SWEtube_density(i,j)==0 %if the quality is bad (0)
               SWEtube_density(i,j-1:j)=nan; %then replace the 0 and the corresponding dneisty value with NaN
            elseif SWEtube_density(i,j)==1 %if the quality is good (1) 
               SWEtube_density(i,j)=nan; %then replace the 1 with a NaN (makes it easier to average later)
            end
        end
        
    end
    
    remove_index = all(isnan(SWEtube_density(:,1:18)),2);
    SWEtube_density(remove_index,:) = [];
    SWEtube_text(remove_index,:) = [];

%Establishing index for various data corresponding to SWEtube_density
    %This was determined by just looking through the imported data
    index_zigzag = [1,3,5:8,10:13,24,26,28:31];
    index_SPSWE = [2,4,9,14,15,25];
    
%Create variables to export for all the SWE value and all snowpit values    
    SWEtube_full = [SWEtube_text, num2cell(SWEtube_density)];
    Snowpit_full = [snowpit_text, num2cell(snowpit_density)];
    depth_full = [depth_text, num2cell(depth_density)];

%% Zigzag Density

%Use index to select only density values from zigzag locations
    zigzagSWE = SWEtube_density(index_zigzag',:);

%Get number of good quality measurements for each location
    count = zeros(size(zigzagSWE,1),1);
    for i = 1:size(zigzagSWE,1) %add together all the 'true' for when a value is present (not nan)
        count(i,1) = sum(~isnan(zigzagSWE(i,1:15)));
    end

%Get mean density and std (kg/m^3)
    % Location, mean density, std, number of obs, elevation
    zigzagSWE = [SWEtube_text(index_zigzag',:), num2cell(nanmean(zigzagSWE(:,1:15),2)),...
        num2cell(nanstd(zigzagSWE(:,1:15),1,2)), num2cell(count), num2cell(SWEtube_density(index_zigzag',21))];

%% Snowpit Density

%Snowpit SWE tube values (based on index)
    snowpittubeSWE = SWEtube_density(index_SPSWE',:);

%Get number of good quality measurements for each location
    count = sum(~isnan(snowpittubeSWE(:,1:15)),2);

%Get mean density and std (kg/m^3)
    % Location, mean density, std, min, max, number of obs
    snowpittubeSWE = [SWEtube_text(index_SPSWE',:), num2cell(nanmean(snowpittubeSWE(:,1:15),2)),...
        num2cell(nanstd(snowpittubeSWE(:,1:15),1,2)), num2cell(min(snowpittubeSWE(:,1:15),[],2)),...
        num2cell(max(snowpittubeSWE(:,1:15),[],2)), num2cell(count)];

%Corresponding snowpit densities
    indexSP = zeros(size(index_SPSWE,2),1);
    for i = 1:size(index_SPSWE,2) %find index where SP label is the same as the snowpit density data labels
        indexSP(i,1) = find(strcmp(snowpittubeSWE(i,1),snowpit_text)); 
    end
    % Location, mean tube density, std, min, max, number of obs, snowpit
    % density, elevation, SP min, SP max
    snowpittubeSWE = [snowpittubeSWE, num2cell(snowpit_density(indexSP,1)), ...
                            num2cell(snowpit_density(indexSP,4:6))];

%% Create structure for relevant data

Density = struct('snowpit',{Snowpit_full},'pitANDtube',{snowpittubeSWE},...
                'tube',{SWEtube_full},'zigzagtube',{zigzagSWE},'SWEdepth',{depth_full});

%% SWE depth correction    

if options.TubeDensity == 2

    % tube = [Density.SWEdepth(4,1), mean(final(4:9,2)), std(final(4:9,2)); Density.SWEdepth(13,1), mean(final(13:19,2)), std(final(13:19,2));...
    %     Density.SWEdepth(29,1), mean(final(29:35,2)), std(final(29:35,2)); Density.SWEdepth(44,1), mean(final(44:50,2)), std(final(44:50,2));...
    %     Density.SWEdepth(54,1), mean(final(54:61,2)), std(final(54:61,2)); Density.SWEdepth(85,1), mean(final(84:89,2)), std(final(84:89,2))];
    % tube(:,2:3) = num2cell(round(cell2mat(tube(:,2:3))));
    % 
    % tubeG = [cellstr('G04'), mean([final(4:9,2);final(13:19,2)]), std([final(4:9,2);final(13:19,2)]);...
    %     cellstr('G02'), mean([final(29:35,2);final(44:50,2)]), std([final(29:35,2);final(44:50,2)]);...
    %     cellstr('G13'), mean([final(54:61,2);final(84:89,2)]), std([final(54:61,2);final(84:89,2)]);...
    %     cellstr('All'), mean([final(4:9,2);final(13:19,2);final(29:35,2);final(44:50,2);final(54:61,2);final(84:89,2)]),...
    %         std([final(4:9,2);final(13:19,2);final(29:35,2);final(44:50,2);final(54:61,2);final(84:89,2)])];
    % tubeG(:,2:3) = num2cell(round(cell2mat(tubeG(:,2:3))));
    % 
    %     x = cell2mat(Density.pitANDtube(:,7)); %Snowpit
    %     y = cell2mat(tube(:,2)); %Detrended SWE tube
    %     errory = [cell2mat(tube(:,3)),cell2mat(tube(:,3))];
    %     %errory = [cell2mat(Density.pitANDtube(:,2))-cell2mat(Density.pitANDtube(:,4)), ...
    %     %    cell2mat(Density.pitANDtube(:,5))-cell2mat(Density.pitANDtube(:,2))]; %min and max tube
    %     errorx = [cell2mat(Density.pitANDtube(:,7))-cell2mat(Density.pitANDtube(:,9)), ...
    %         cell2mat(Density.pitANDtube(:,10))-cell2mat(Density.pitANDtube(:,7))]; %min and max SP
    % errorbarxy(x,y,errorx(:,1),errorx(:,2),errory(:,1),errory(:,2),'Color','k','LineStyle','none','Marker','o',...
    %    'MarkerFaceColor','k','LineWidth',1,'MarkerSize',5); hold on
    %     
    %     P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    %     LM = fitlm(x,y);
    % plot(x,yfit,'r')
    %     xlabel('Snowpit density (kg m^{-3})')
    %     ylabel('SWE tube density (kg m^{-3})')
    %     dim = [0.15,0.5,0.11,0.11];
    %     str = {strcat('y=',num2str(round(P(1),2)),'*x+',num2str(round(P(2),2))), ...
    %         strcat('R^2=',num2str(round(LM.Rsquared.Ordinary,2)))} ;
    %     annotation('textbox',dim,'String', str,'FitBoxToText','on')
    %     text(x+2, y+2, Density.pitANDtube(:,1))
    %      axis([220 400 220 400])
    %     axis equal 


        depth = cell2mat(Density.SWEdepth(:,3));
        den = cell2mat(Density.SWEdepth(:,4));
    C = [depth, den];
    [~, I] = sortrows(C);
    detrendC = [I, detrend(C(:,2)) + mean(C(:,2))];
    final = sortrows(detrendC);

    % plot(detrendC(:,2),detrendC(:,1),'.','markers',13)
    %     xlabel('Snowpit density (kg m^{-3})')
    %     ylabel('SWE tube density (kg m^{-3})')

    final2 = [Density.SWEdepth(:,1), num2cell(final(:,2))];
    a = categorical(final2(:,1));
    a1 = [a(1:end-1),a(2:end)];
    ind = [find(a1(:,1)~=a1(:,2));length(a)];
    index = [1;ind(1:end-1)+1];
    index = [index, ind];
    M = cell(10,6);
    for i = 1:length(index)
        M(i,1:6) =  [final2(index(i,1),1), num2cell(mean(cell2mat(final2(index(i,1):index(i,2),2)))),...
            num2cell(std(cell2mat(final2(index(i,1):index(i,2),2)))),...
            num2cell(min(cell2mat(final2(index(i,1):index(i,2),2)))),...
            num2cell(max(cell2mat(final2(index(i,1):index(i,2),2)))),...
            (index(i,2)-index(i,1)+1)];
    end

    Density.tube = [M, Density.tube(:,20:22)];
    Density = rmfield(Density, 'SWEdepth');
    Density.pitANDtube(:,2:6) = M(index_SPSWE,2:6);
    Density.zigzagtube(:,2) = M(index_zigzag,2);

end
%% Clear variables
clear SWEtube* snowpit* zigzagSWE Snowpit_full depth*
clear P LM x y index* i j count yfit str dim filename error
clear a* C den detrendC final* I ind M remove_index