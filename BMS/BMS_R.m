function [ BMSbest ] = BMS_R( SWE, topo_sampled_ns )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

run /Users/Alexandra/Documents/SFU/Data/SnowDepth/OPTIONS.m
options.ZZ = 2; %exclude zigzags
run /Users/Alexandra/Documents/SFU/Data/SnowDepth/MAIN

%% Make Cal Val data sets
for g = 1:3
    glacier = char(options.glacier(g));
    y = SWE(g).swe;

 %Choose number of runs
runs = 1000;        

 %Cross validation random number matrix
[~, cal_ind] = sort(rand(runs,length(y)),2); %create matrix of random numbers
cal_ind = cal_ind(:,1:floor(length(y)*3/4)); %Choose 3/4 for the calibration component


% Run BMS code in R with cross validation

for i = 1:runs                              %for number of runs
        cal_ind_temp    = cal_ind(i,:);         %get random numbers for choosing obs
        val_ind         = setdiff(1:length(y),cal_ind_temp); %get random numbers for validating

        sweG = SWE(g).swe(cal_ind_temp);

        heads   = fieldnames(topo_sampled_ns.G4);
        Vtopo.(glacier) = []; 
        for j = 1:length(heads)
            param             = char(heads(j));
            topoG.(param)  = topo_sampled_ns.(glacier).(param)(cal_ind_temp);
            Vtopo.(glacier) = [Vtopo.(glacier), topo_sampled_ns.(glacier).(param)(val_ind)];
        end

        save mat2R.mat sweG topoG
        % Run BMS code in R
        %!R CMD BATCH BMS_matlab.R
        !/usr/local/bin/R CMD BATCH BMS_matlab.R
        
        % Make table of coeffs
            load R2mat.mat

            tableCol = {'Coefficient','SD', 'PIP','PSP'};
            tableRow = {'aspect','elevation','northness','profileCurv','slope','tangentCurve','Sx','centreD','intercept'};

            BMSg.uniform(i,g) = {table(Gcoeffs.Post_Mean,  Gcoeffs.Post_SD, Gcoeffs.PIP, Gcoeffs.Cond_Pos_Sign,...
                                'RowNames', tableRow, 'VariableNames',tableCol)};
            BMSg.fixed(i,g) = {table(Gcoeffs.Post_Mean_1, Gcoeffs.Post_SD_1, Gcoeffs.PIP_1, Gcoeffs.Cond_Pos_Sign_1,...
                                'RowNames', tableRow, 'VariableNames',tableCol)};
            BMSg.random(i,g) = {table(Gcoeffs.Post_Mean_2,Gcoeffs.Post_SD_2,  Gcoeffs.PIP_2, Gcoeffs.Cond_Pos_Sign_2,...
                                'RowNames', tableRow, 'VariableNames',tableCol)};                       
            BMSg.variable(i,g) = {table(Gcoeffs.Post_Mean_3, Gcoeffs.Post_SD_3, Gcoeffs.PIP_3, Gcoeffs.Cond_Pos_Sign_3,...
                                'RowNames', tableRow, 'VariableNames',tableCol)};        
            %BMSg.mcmc = table(Gcoeffs.Post_Mean_4, Gcoeffs.Post_SD_4, Gcoeffs.PIP_4, Gcoeffs.Cond_Pos_Sign_4,...
            %                    'RowNames', tableRow, 'VariableNames',tableCol);        
                         
            models = fieldnames(BMSg);
            
           % validation
            for m = 1:length(models)
                mod = char(models(m));
                
                y_regress = sum(Vtopo.(glacier).*repmat(BMSg.uniform{i,g}{1:end-1,1}',...
                                    length(Vtopo.(glacier)),1),2) + BMSg.uniform{i,g}{end,1};      
                rmse.(mod).(glacier)(i,1) = sqrt(sum((SWE(g).swe(val_ind)-...
                                            y_regress).^2)/numel(y_regress));
            end
                        
            
end

    %min rmse coefficients
    BMSbest.(glacier) = table();
    for m = 1:length(models)
        mod = char(models(m));
        
        rmse_min.(glacier)    = rmse.(mod).(glacier)==min(rmse.(mod).(glacier));
        bestCoefs   = [BMSg.(mod){rmse_min.(glacier),g}(:,1); ...
            table(min(rmse.(mod).(glacier)),'VariableNames',{'Coefficient'},'RowNames',{'rmse'})]; 
        bestCoefs.Properties.VariableNames = {mod};
        
        BMSbest.(glacier) = [BMSbest.(glacier), bestCoefs];
    end
end    
   
end

