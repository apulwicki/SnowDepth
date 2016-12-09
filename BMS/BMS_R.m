function [ BMSbest, residuals ] = BMS_R( SWE, topoSampled )
% Bayesian model averaging with cross validation 
%       This script uses a BMA R package to compute coefficients for the 
%       linear regression with cross validation for some number
%       of runs (typically 1000). The model with the lowest RMSE is chosen.
%       The percent variance explained by each coefficient and
%       the residuals of the predicted data are also calculated.
%
%       Inputs:         SWE values (SWE)
%                       Topographic param structure (topoSampled)
%       Outputs:        BMA coefficients and % var explained (BMSbest)
%                       Residuals of predicted data (residuals)

%       Alexandra Pulwicki  Created: December 2016
%                           Updated: December 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for g = 1:length(SWE) %For all glaciers
%% Make Cal Val data sets
 %Initialize
    GG = {'G4','G2','G13'};
    glacier = char(GG(g));
    y = SWE(g).swe;

 %Choose number of runs
runs = 10;        

 %Cross validation random number matrix
[~, cal_ind] = sort(rand(runs,length(y)),2);      %create matrix of random numbers
cal_ind      = cal_ind(:,1:floor(length(y)*3/4)); %Choose 3/4 for the calibration component

 % Run BMS code in R with cross validation
for i = 1:runs                                  %for number of runs
        cal_ind_temp    = cal_ind(i,:);         %get random numbers for choosing obs
        val_ind         = setdiff(1:length(y),cal_ind_temp); %get random numbers for validating

        sweG = SWE(g).swe(cal_ind_temp);        %select calibration swe values

         %Select calibration and validation topographic data
            heads           = fieldnames(topoSampled.G4);   %get parameter names
            Vtopo.(glacier) = [];                           %initalize validation matrix
        for j = 1:length(heads)
            param           = char(heads(j));           
            topoG.(param)   = topoSampled.(glacier).(param)(cal_ind_temp);  %calibration topo
            Vtopo.(glacier) = [Vtopo.(glacier), topoSampled.(glacier).(param)(val_ind)]; %validation topo
        end
         %Save the swe and topo variables 
        save mat2R.mat sweG topoG
        
         %Run BMS code in R though the terminal
        %!R CMD BATCH BMS_matlab.R
        !/usr/local/bin/R CMD BATCH BMS_matlab.R
        
         %Load data from R
        load R2mat.mat                                  
         
         %Make table of coefficients for each run
            tableCol = {'Coefficient','SD', 'PIP','PSP'};   %variable names
            tableRow = [heads; {'Intercept'}];              %parameters names as row names
        BMSg(i,g) = {table(Gcoeffs.Post_Mean,  Gcoeffs.Post_SD, Gcoeffs.PIP, Gcoeffs.Cond_Pos_Sign,...
                             'RowNames', tableRow, 'VariableNames',tableCol)}; %create table with coefficient values for each run
         %RMSE
        y_regress = sum(Vtopo.(glacier).*repmat(BMSg{i,g}{1:end-1,1}',...
                       length(Vtopo.(glacier)),1),2) + BMSg{i,g}{end,1};  %modelled values   
        rmse.(glacier)(i,1) = sqrt(sum((SWE(g).swe(val_ind)-...
                        y_regress).^2)/numel(y_regress)); %RMSE between modelled and observed

        %BMS coeffs with lowest RMSE
        rmse_min.(glacier)    = rmse.(glacier)==min(rmse.(glacier));  %min RMSE value
        BMSbest.(glacier)             = [BMSg{rmse_min.(glacier),g}(:,1); ... %coeffs for that run
                table(min(rmse.(glacier)),'VariableNames',{'Coefficient'},'RowNames',{'rmse'})];         
end    
   
%% Calculate % variance explained by each variable

beta    = BMSbest.(glacier).Properties.RowNames(1:end-2);   %names of params
SSt     = sumsqr(y-mean(y));                %total sum of squares

Pvar = table(zeros(length(beta),1),'RowNames',beta);        %initalize
Pvar.Properties.VariableNames = {'PercentVarExplained'};
for i = 1:length(beta)                      %only coeffs, no intercept
    rowname         = char(beta(i));        %coeff name
    Xfit            = BMSbest.(glacier){end-1,1} + ...
                        BMSbest.(glacier){i,1}*topoSampled.(glacier).(rowname); %topo params times their coeffs
    SSr             = sumsqr(Xfit-mean(y)); %residual sum of squares
    Pvar{i,1}       = SSr/SSt*100;          %percent var exmaplined
end

Pvar = [Pvar; table([0; 0], 'RowNames',BMSbest.(glacier).Properties.RowNames(end-1:end),...
                            'VariableNames',{'PercentVarExplained'})];
BMSbest.(glacier) = [BMSbest.(glacier), Pvar];              %add to final table

%% Residuals
 %Predict all data at sampling locations
X1 = struct2table(topoSampled.(glacier));   %topo data
y_regress = sum(X1{:,:}.*...                %modelled swe
    repmat(BMSbest.(glacier){1:end-2,1}',height(X1),1),2) + BMSbest.(glacier){end-1,1};      %predict validation data

 %Residual
residuals.(glacier) = y-y_regress;
end

