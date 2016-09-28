function [mlr_best, rmse_best] = MLRcalval(y, X)

n = size(X,2);
Xones = [ones(length(X),1), X];
c = [ones(2^n,1), (dec2bin(0:(2^n)-1)=='1')];      c = c(2:end,:);
mlr_best = zeros(length(c),1);    rmse_best = zeros(length(c),1);

for j = 1:length(c)
    mult = repmat(c(j,:),length(X),1);
    X1 = Xones.*mult; X1(:,~any(X1,1)) = [];  

    mDim = size(X1,2);
    runs = 1000;
    rmse = zeros(runs,1); mlr = zeros(runs,mDim);

    for i = 1:runs
        cal_ind = randperm(length(y),floor(length(y)/2));
        val_ind = setdiff(1:length(y),cal_ind);

        mlr(i,:) = regress(y(cal_ind,1),X1(cal_ind,:))';

        y_regress = sum(repmat(mlr(i,:),length(val_ind),1).*X1(val_ind,:),2);
        rmse(i,1) = sqrt(sum((y(val_ind,1)-y_regress).^2)/numel(y_regress));

    end

    %plot(rmse,'.')
    
    col = find(c(j,:));
    for i = 1:size(col,2)
        mlr_best(j,col(1,i)) = mlr(rmse==min(rmse),:);
    end
    rmse_best(j,1)  = min(rmse);
end
end