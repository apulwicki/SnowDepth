function [mlr_final, rmse_final, lm] = MLRcalval(y, X)

%Convert from structure to matrix
params = fieldnames(X); 
M = X.(char(params(1)));
for j = 2:length(params)
    field = char(params(j));
    M = [M, X.(field)];
end

n = size(M,2);
Xones = [ones(length(M),1), M];
c = [ones(2^n,1), (dec2bin(0:(2^n)-1)=='1')];      c = c(2:end,:);
mlr_best = zeros(length(c),1);    rmse_best = zeros(length(c),1); 

runs = 1000;        
[~, cal_ind] = sort(rand(runs,length(y)),2);
cal_ind = cal_ind(:,1:floor(length(y)*3/4));

for j = 1:length(c)
    mult = repmat(c(j,:),length(M),1);
    X1 = Xones.*mult; X1(:,~any(X1,1)) = [];  

    mDim = size(X1,2);
    rmse = zeros(runs,1); mlr = zeros(runs,mDim);

    for i = 1:runs
        cal_ind_temp = cal_ind(i,:);
        val_ind = setdiff(1:length(y),cal_ind_temp);

        mlr(i,:) = regress(y(cal_ind_temp,1),X1(cal_ind_temp,:))';

        y_regress = sum(repmat(mlr(i,:),length(val_ind),1).*X1(val_ind,:),2);
        rmse(i,1) = sqrt(sum((y(val_ind,1)-y_regress).^2)/numel(y_regress));

    end

    %plot(rmse,'.')
    
    col = find(c(j,:));
    for i = 1:size(col,2)
        mlr_best(j,col(1,i)) = mlr(rmse==min(rmse),i);  %chose best coefficients  
    end
    rmse_best(j,1)  = min(rmse); %lowest rmse value
end

%Chose best combination of variables and output final MLR coefficients
best = rmse_best==min(rmse_best); 

%Do a linear regression and get the stats (esp. p value) for the
%coefficents
c_best = c(best,2:end); c_best = repmat(c_best,length(M),1); 
T = M.*c_best;  T = T(:,any(T));
lm = fitlm(T,y,'linear'); 

%Caluclate % varience explained by each variable
an = anova(lm); SumSq = table2array(an(:,1));
Pvar = [0; SumSq(1:end-1)/sum(SumSq)*100]; Pvar = table(Pvar);

coeffs = [lm.Coefficients(:,1),lm.Coefficients(:,4), Pvar];
rmse_final = rmse_best(best,1);

mlr_final = table(zeros(8,1),zeros(8,1),zeros(8,1),...
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




