% Importing Snow Density Data
%       This script imports snow density data from the summary of density
%       data file. It also differentes the SWE tube density values and the
%       snowpit density values.

%       Inputs:         Summary of density data ('summary_densitydata.xlsx')
%       Outputs:        Density data structure with separated values for
%                           just snowpit, snowpit and tube, just tube, and
%                           tube in zigzag locations. 

%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing density values

%Get data from file
    [snowpit_density, snowpit_text, snowpit_raw] = xlsread('summary_densitydata.xlsx','snowpit','A2:E11');
    [SWEtube_density, SWEtube_text, SWEtube_raw] = xlsread('summary_densitydata.xlsx','SWEtube','A2:V34');

%Establishing index for various data corresponding to SWEtube_density
    %This was determined by just looking through the imported data
    index_zigzag = [1,3,5,8,10,13,26,28,30,32];
    index_SPSWE = [2,4,9,14,15,27];

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

%Create variables to export for all the SWE value and all snowpit values    
    SWEtube_full = [SWEtube_text, num2cell(SWEtube_density)];
    Snowpit_full = [snowpit_text, num2cell(snowpit_density)];
    
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
    count = zeros(size(snowpittubeSWE,1),1);
    for i = 1:size(snowpittubeSWE,1)%add together all the 'true' for when a value is present (not nan)
        count(i,1) = sum(~isnan(snowpittubeSWE(i,1:15)));
    end

%Get mean density and std (kg/m^3)
    % Location, mean density, std, number of obs
    snowpittubeSWE = [SWEtube_text(index_SPSWE',:), num2cell(nanmean(snowpittubeSWE(:,1:15),2)),...
        num2cell(nanstd(snowpittubeSWE(:,1:15),1,2)), num2cell(count)];

%Corresponding snowpit densities
    indexSP = zeros(size(index_SPSWE,2),1);
    for i = 1:size(index_SPSWE,2) %find index where SP label is the same as the snowpit density data labels
        indexSP(i,1) = find(strcmp(snowpittubeSWE(i,1),snowpit_text)); 
    end
    % Location, mean density, std, number of obs, snowpit density, elevation
    snowpittubeSWE = [snowpittubeSWE, num2cell(snowpit_density(indexSP,1)), ...
                            num2cell(snowpit_density(indexSP,4))];

%% Create structure for relevant data

Density = struct('snowpit',{Snowpit_full},'pitANDtube',{snowpittubeSWE},...
                'tube',{SWEtube_full},'zigzagtube',{zigzagSWE});
                        
%% Clear variables
clear SWEtube* snowpit* zigzagSWE Snowpit_full
clear P LM x y index* i j count yfit str dim filename error
