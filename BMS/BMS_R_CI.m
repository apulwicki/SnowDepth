function [ BMSbest, residuals ] = BMS_R_CI( swe, topoSampled )
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
global options

for g = 1:3 %For all glaciers
%% Make Cal Val data sets

 %Initialize
    glacier = options.glacier{g}; display(['glacier = ', glacier]);
    y = swe.(glacier)(:,1);

 %Choose number of runs
runs = 100;        

 %Cross validation random number matrix
[~, cal_ind] = sort(rand(runs,length(y)),2);      %create matrix of random numbers
cal_ind      = cal_ind(:,1:floor(length(y)*2/3)); %Choose 2/3 for the calibration component

 % Run BMS code in R with cross validation
for i = 1:runs                                  %for number of runs
        cal_ind_temp    = cal_ind(i,:);         %get random numbers for choosing obs
        val_ind         = setdiff(1:length(y),cal_ind_temp); %get random numbers for validating

        sweG = y(cal_ind_temp);        %select calibration swe values

         %Select calibration and validation topographic data
            heads           = fieldnames(topoSampled.G4);   %get parameter names
            Vtopo.(glacier) = [];                           %initalize validation matrix
        for j = 1:length(heads)
            param           = char(heads(j));           
            topoG.(param)   = topoSampled.(glacier).(param)(cal_ind_temp);  %calibration topo
            Vtopo.(glacier) = [Vtopo.(glacier), topoSampled.(glacier).(param)(val_ind)]; %validation topo
        end
         %Save the swe and topo variables 
        save mat2R_CI.mat sweG topoG
        
         %Run BMS code in R though the terminal (top SFU, bottom Mac)
        !R --no-save CMD BATCH BMS_matlab_CI.R
%         !/usr/local/bin/R CMD BATCH BMS_matlab.R    
        
        %Load data from R
        try
            load R2mat_CI.mat 
            delete R2mat_CI.mat
        catch
            !R --no-save CMD BATCH BMS_matlab_CI.R
            load R2mat_CI.mat 
            delete R2mat_CI.mat
        end
         
         %Make table of coefficients for each run
            tableCol = {'Coefficient','SD', 'PIP','PSP'};   %variable names
            tableRow = [heads; {'Intercept'}];              %parameters names as row names
        BMSg(i,g) = {table(Gcoeffs.Post_Mean,  Gcoeffs.Post_SD, Gcoeffs.PIP, Gcoeffs.Cond_Pos_Sign,...
                             'RowNames', tableRow, 'VariableNames',tableCol)}; %create table with coefficient values for each run
         %RMSE
        y_regress = sum(Vtopo.(glacier).*repmat(BMSg{i,g}{1:end-1,1}',...
                       size(Vtopo.(glacier),1),1),2) + BMSg{i,g}{end,1};  %modelled values   
        rmse.(glacier)(i,1) = sqrt(sum((y(val_ind)-...
                        y_regress).^2)/numel(y_regress)); %RMSE between modelled and observed

        %BMS coeffs with lowest RMSE
        rmse_min.(glacier)    = find(rmse.(glacier)==min(rmse.(glacier)));  %min RMSE value
        BMSbest.(glacier)     = [BMSg{rmse_min.(glacier)(1,1),g}(:,1); ... %coeffs for that run
                table(min(rmse.(glacier)),'VariableNames',{'Coefficient'},'RowNames',{'rmse'})];         
end    
   
%% Calculate % variance explained by each variable

beta    = BMSbest.(glacier).Properties.RowNames(1:end-2);   %names of params

%--------Semi-partial (Part) correlation squared
    %Code adapted from "ppcor: An R Package for a Fast Calculation to  
    %       Semi-partial Correlation Coefficients" Kim 2015
semiR = table(zeros(length(beta),1),'RowNames',beta);    %initalize
semiR.Properties.VariableNames = {'SemiR2'};

    M = struct2table(topoSampled.(glacier));
try
cx = cov([y,M{:,:}]);   %Covariance matrix of data
dx = inv(cx);           %Inverse covariance matrix
pc = -corrcov(dx);      %Convert correlations to covariance
    n = size(pc,1);     pc(1:(n+1):end) = 1; %Set diagonal elements to 1

kk  = pc./repmat(sqrt(diag(cx)),1,n)./...   %Semi-partial correlation (Eq. 2.6) 
    sqrt(abs(repmat(diag(dx),1,n)-((dx.^2).'./repmat(diag(dx),1,n)).'));
kk(1:(n+1):end) = 1;    %Set diagonal elements to 1
kk = kk.^2;             %Square correlations

semiR{:,1} = kk(1,2:end)';  %Assign to table
    %Alternative method = calculate residuals of var of interest with other
    %vars and then correlate residuals with y {res = fitlm([deg',disp'],BC); 
    %corr(hl',res.Residuals.Raw)}
catch
   semiR{:,1} =  NaN([length(beta),1]);
   disp('Error handled')
end
    
%--------Univariate R-squared
uniR = table(zeros(length(beta),1),'RowNames',beta);    %initalize
uniR.Properties.VariableNames = {'UnivarR2'};

uniR{:,1}   = corr(M{:,:},y).^2;      %Squared raw correlation between
                                           %regressors and y data

 %Add to final table
Pvar = [semiR, uniR];        
Pvar = [Pvar; table([0; 0], [0; 0], 'RowNames',BMSbest.(glacier).Properties.RowNames(end-1:end),...
                            'VariableNames',Pvar.Properties.VariableNames)];
BMSbest.(glacier) = [BMSbest.(glacier), Pvar];              %add to final table
BMSbest.(glacier) = [BMSbest.(glacier); table(0,0,0, 'RowNames',{'R2'},...
                            'VariableNames',BMSbest.(glacier).Properties.VariableNames)];

%% Residuals
 %Predict all data at sampling locations
y_regress = sum(M{:,:}.*...                %modelled swe
    repmat(BMSbest.(glacier){1:end-3,1}',height(M),1),2) + BMSbest.(glacier){end-2,1};      %predict validation data

 %Residual
residuals.(glacier) = y-y_regress;

 %R2 of fit
BMSbest.(glacier){end,2} = corr(y_regress,y)^2;                        

end

