function [mlr_final, rmse_final, lm] = MLRcalval(y, X)

%Convert from structure to table
M = struct2table(X);

n = size(M,2);
c = logical(dec2bin(0:(2^n)-1)=='1');      c = c(2:end,:);
mlr_best = zeros(length(c),1);    rmse_best = zeros(length(c),1); 

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
    coeffs_w{j,1} = mlr_best{1,1}.Coefficients(:,1);     
    coeffs_w{j,1}{:,1} = coeffs_w{j,1}{:,1}*AICweight(1,1);
end


%Do a linear regression and get the stats (esp. p value) for the
%coefficents
c_best = c(best,2:end); c_best = repmat(c_best,length(M),1); 
T = M.*c_best;  T = T(:,any(T));
lm = fitlm(T,y,'linear'); 

%Caluclate % varience explained by each variable
an = anova(lm); SumSq = table2array(an(:,1));
Pvar = [0; SumSq(1:end-1)/sum(SumSq)*100]; Pvar = table(Pvar);

%AIC values
aic = lm.ModelCriterion.AIC; %AIC value

coeffs = [lm.Coefficients(:,1),lm.Coefficients(:,4), Pvar, aic];

    rmse_final = rmse_best(best,1);

mlr_final = table(zeros(length(params)+1,1),zeros(length(params)+1,1),zeros(length(params)+1,1),...
    'RowNames',[{'intercept'};params],...
    'VariableNames',{'coefficient','pvalue','PercentVarExplained'});
mlr_final(1,1:2) = coeffs(1,1:2);  

row = 2; next = 2;
for i = 1:length(c_best(1,:))
    if c_best(1,i)~=0
       mlr_final(row,1:3) = coeffs(next,1:3);
       next = next+1;
    end
    row = row+1;
end        
    
end




