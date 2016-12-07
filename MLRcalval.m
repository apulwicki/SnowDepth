function [coeffs_final, residuals] = MLRcalval(y, X)
% MLR with cross validation and linear combos of params
%       This script calculates a MLR with cross validation for some number
%       of runs (typically 1000) for all possible linear combination of
%       parameters (models). The models are then weighted using BIC values
%       and the output is a single coefficient for each topographic
%       parameter. The percent variance explained by each coefficient and
%       the residuals of the predicted data are also calculated
%
%       Inputs:         SWE values (y)
%                       Topographic param structure (X)
%       Outputs:        MLR coefficients and % var explained (coeffs_final)
%                       Residuals of predicted data (residuals)

%       Alexandra Pulwicki  Created: September 2016
%                           Updated: November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Initializing function

 %Convert from structure to table
M = struct2table(X);
n = size(M,2);

 %Make a logical matrix with all linear combinations for choosing
 %parameters
c = logical(dec2bin(0:(2^n)-1)=='1');      c = c(2:end,:);

 %Choose number of runs
runs = 1000;        

 %Cross validation random number matrix
[~, cal_ind] = sort(rand(runs,length(y)),2); %create matrix of random numbers
cal_ind = cal_ind(:,1:floor(length(y)*3/4)); %Choose 3/4 for the calibration component

 %Initialize matrices
mlr_best = cell(length(c),1);    rmse_best = zeros(length(c),1);    BIC_best = rmse_best;

%% MLR with cross validation

for j = 1:length(c)                             %all linear combos of params
    X1 = M(:,c(j,:));                           %create matrix with one combo of params
    rmse = zeros(runs,1);   mlr = cell(runs,1); %initialize output matrices

    for i = 1:runs                              %for number of runs
        cal_ind_temp    = cal_ind(i,:);         %get random numbers for choosing obs
        val_ind         = setdiff(1:length(y),cal_ind_temp); %get random numbers for validating
        
        mlr{i,1}        = regress(y(cal_ind_temp,1), ...     %MLR between calibration data and chosen topo params
            [ones(length(cal_ind_temp),1), X1{cal_ind_temp,:}])'; %(need to include a row of ones for getting intercept)
        
        y_regress   = sum(X1{val_ind,:}.*mlr{i,1}(2:end),2) + mlr{i,1}(1);      %predict validation data
        rmse(i,1)   = sqrt(sum((y(val_ind,1)-y_regress).^2)/numel(y_regress));  %get the RMSE between observed and predicted values

    end

    %plot(rmse,'.')    %hist(rmse)
    %[h, p]=chi2gof(rmse);    display(['j = ',num2str(j),' h = ',num2str(h), ' p = ',num2str(p)])
    
    rmse_best(j,1)  = min(rmse);                        %lowest rmse value

     %Do fitlm for best set of calibration data
    min_rmse        = rmse==min(rmse);                  %find best set
    swe_obs         = table(y(cal_ind(min_rmse,:)), 'VariableNames',{'swe'}); %create table of input data
    topo_obs        = M(cal_ind(min_rmse,:),c(j,:));    %create table of topo params
    mlr_best{j,1}   = fitlm([topo_obs, swe_obs]);       %do thorough MLR using fitlm
    BIC_best(j,1)   = mlr_best{j,1}.ModelCriterion.BIC; %BIC of best MLR run

end

%% Weighting models 

 %Calculate weight for each model based on BIC value
BICweight = exp(-(BIC_best-min(BIC_best))/2);   %exponential weighting compared to best (min) BIC
BICweight = BICweight/sum(BICweight);           %normlize weights

 %Get the coefficients for each param and weight them 
coeffs_w{j,1} = cell(length(mlr_best),1);                       %initialize
for j = 1:length(mlr_best)
    coeffs_w{j,1} = mlr_best{j,1}.Coefficients(:,1);            %get coeffs
    coeffs_w{j,1}.Properties.VariableNames = {['C',num2str(j)]};%name column with combo number
    coeffs_w{j,1}{:,1} = coeffs_w{j,1}{:,1}*BICweight(j,1);     %weights the coefficients
end

 %Fill in missing coefficients (the zeors ones)
all = coeffs_w{end,1}.Properties.RowNames;                      %get param names 
coeffs_full = table();                                          %initialize
for j = 1:length(mlr_best)
    missing = all(~ismember(all,coeffs_w{j,1}.Properties.RowNames)); %determine which ones are zeros
    T = table(zeros(length(missing),1),'RowNames',missing,...   %create a table with full coefficients
        'VariableNames',coeffs_w{j,1}.Properties.VariableNames);
    coeffs_w{j,1} = [coeffs_w{j,1};T];                          

    coeffs_full = [coeffs_full,coeffs_w{j,1}];                  %append to the table
end

 %Sum over all weighted coefficients and return nice table with final
 %coeffs
for i = 1:height(coeffs_full)
    coeffs_final(i,1) = table(sum(coeffs_full{i,:}));
end
coeffs_final.Properties.RowNames = coeffs_full.Properties.RowNames;
coeffs_final.Properties.RowNames(1,1) = {'Intercept'};
coeffs_final.Properties.VariableNames = {'Coefficient'};


%% Caluclate % variance explained by each variable

beta    = coeffs_final.Properties.RowNames; %names of params
SSt     = sumsqr(y-mean(y));                %total sum of squares

Pvar = table(zeros(length(beta),1),'RowNames',beta);    %initalize
Pvar.Properties.VariableNames = {'PercentVarExplained'};
for i = 2:length(beta)                      %only coeffs, no intercept
    rowname         = char(beta(i));        %coeff name
    Xfit            = coeffs_final{1,1} + coeffs_final{i,1}*X.(rowname); %topo params times their coeffs
    SSr             = sumsqr(Xfit-mean(y)); %residual sum of squares
    Pvar{i,1}       = SSr/SSt*100;          %percent var exmaplined
end

coeffs_final = [coeffs_final, Pvar];        %add to final table

%% Find rmse of coefficients

 %Predict all data at sampling locations
y_regress   = sum(X1{:,:}.*repmat(coeffs_final{:,1}(2:end)',height(X1),1),2) + coeffs_final{1,1};      %predict validation data
rmse_final  = sqrt(sum((y-y_regress).^2)/numel(y_regress));  %get the RMSE between observed and predicted values

coeffs_final = [coeffs_final(2:end,:); coeffs_final(1,:)];
coeffs_final = [coeffs_final; table(rmse_final, 0, ...
                'VariableName',{'Coefficient', 'PercentVarExplained'},'RowNames',{'rmse'})];

%% Residuals

 %Get residual
residuals = y-y_regress;

end




