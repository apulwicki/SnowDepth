function [coeffs_final, residuals] = MLRcalval(y, X)

%Convert from structure to table
M = struct2table(X);

n = size(M,2);
c = logical(dec2bin(0:(2^n)-1)=='1');      c = c(2:end,:);
mlr_best = cell(length(c),1);    rmse_best = zeros(length(c),1); 

runs = 1000;        
[~, cal_ind] = sort(rand(runs,length(y)),2);
cal_ind = cal_ind(:,1:floor(length(y)*3/4));

for j = 1:length(c)
    X1 = M(:,c(j,:));  

    rmse = zeros(runs,1); mlr = cell(runs,1);

    for i = 1:runs
        cal_ind_temp = cal_ind(i,:);
        val_ind = setdiff(1:length(y),cal_ind_temp);
        swe_obs = table(y(cal_ind_temp,1), 'VariableNames',{'swe'});
        X2 = [X1(cal_ind_temp,:), swe_obs];
        
        mlr{i,1} = fitlm(X2);
        
        y_regress = predict(mlr{i,1},X1(val_ind,:));
        rmse(i,1) = sqrt(sum((y(val_ind,1)-y_regress).^2)/numel(y_regress));

    end

    %plot(rmse,'.')
    %hist(rmse)
    %[h, p]=chi2gof(rmse);    display(['j = ',num2str(j),' h = ',num2str(h), ' p = ',num2str(p)])

    mlr_best{j,1}   = mlr{rmse==min(rmse)};  %choose best coefficients  
    rmse_best(j,1)  = min(rmse); %lowest rmse value
    AIC_best(j,1)   = mlr_best{j,1}.ModelCriterion.AIC;
end

AICweight = exp(-(AIC_best-min(AIC_best))/2); AICweight = AICweight/sum(AICweight);

for j = 1:length(mlr_best)
    coeffs_w{j,1} = mlr_best{j,1}.Coefficients(:,1); 
    coeffs_w{j,1}.Properties.VariableNames = {['C',num2str(j)]};
    coeffs_w{j,1}{:,1} = coeffs_w{j,1}{:,1}*AICweight(j,1);
end

all = coeffs_w{end,1}.Properties.RowNames;
coeffs_full = table();      
for j = 1:length(mlr_best)
    missing = all(~ismember(all,coeffs_w{j,1}.Properties.RowNames));
    T = table(zeros(length(missing),1),'RowNames',missing,'VariableNames',coeffs_w{j,1}.Properties.VariableNames);
    coeffs_w{j,1} = [coeffs_w{j,1};T];

    coeffs_full = [coeffs_full,coeffs_w{j,1}];
end

for i = 1:height(coeffs_full)
    coeffs_final(i,1) = table(sum(coeffs_full{i,:}));
end
coeffs_final.Properties.RowNames = coeffs_full.Properties.RowNames;
coeffs_final.Properties.RowNames(1,1) = {'Intercept'};
coeffs_final.Properties.VariableNames = {'Coefficient'};


%% Caluclate % variance explained by each variable

beta    = coeffs_final.Properties.RowNames;
SSt     = sumsqr(y-mean(y));

Pvar = table(zeros(length(beta),1),'RowNames',beta);
Pvar.Properties.VariableNames = {'PercentVarExplained'};
for i = 2:length(beta)
    rowname         = char(beta(i));
    Xfit            = coeffs_final{1,1} + coeffs_final{i,1}*X.(rowname);
    SSr             = sumsqr(Xfit-mean(y));
    Pvar{i,1}       = SSr/SSt*100;
end

coeffs_final = [coeffs_final, Pvar];

%% Residuals

rows = coeffs_final.Properties.RowNames;
for i = 1:length(rows)
    param = char(rows(i));
    B.(param) = coeffs_final{i,1};
end

Yfit = repmat(B.Intercept, length(X.(param)),1);
for i = 2:length(rows)
    param = char(rows(i));
    Yfit =  B.(param)*X.(param) + Yfit;
end

residuals = y-Yfit;

end




