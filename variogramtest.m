%test variogram

    %pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all','G04','all','LH',1,'fat');

data = [nanmean(z(5).depth(:,1:4),2), z(5).depth(:,6:7)];

%% Calculate distance between points (lag)

%find the variance for each pair of points
variotest = zeros(size(data,1)*size(data,1),1);
for i = 1:size(data,1)
    z = data(i,1); x = data(i,2); y= data(i,3);
    for j = 1:size(data,1)
        variotest(i*j,1) = EuclideanDistance(x,y,data(j,2),data(j,3)); 
        variotest(i*j,2) = (z-data(j,1))^2;
    end
end

lag = 10;
maxlag = round(max(variotest(:,1)),-1);
intervals = [lag:lag:maxlag]';

variodata = zeros(maxlag/lag,2);
for n = lag:lag:maxlag+lag
    index = intersect(find(variotest(:,1)>n), find(variotest(:,1)<n+1));
    variodata(n/lag,1:2) = [n, sum(variotest(index,2))/(2*length(index))];
end
n = [0, lag];
index = intersect(find(variotest(:,1)>n(1,1)), find(variotest(:,1)<n(1,2)));
variodata = [[n(1,1), sum(variotest(index,2))/(2*length(index))]; variodata];


scatter(variodata(:,1),variodata(:,2))

clear x y z n j i index variotest
%% Variogram function

z = pulldata(SD,'all','G04','all','LH',1,'fat');
    x = z(5).depth(:,6); x2 = nanmax(x)-x; %convert easting to distance in m
    y = z(5).depth(:,7); y2 = nanmax(y)-y; %convert easting to distance in m
    z = nanmean(z(5).depth(:,1:4),2);

[d, ~]= variogram([x2 y2],z,'plotit',false,'nrbins',200);

figure 
        h=d.distance;
        gammaexp = d.val;
        numobs = d.num;
        a0 = 500; % initial value: range 
        c0 = 500; % initial value: sill     
        [a,b,n] = variogramfit(h,gammaexp,a0,c0,numobs,...
                               'solver','fminsearchbnd',...
                               'nugget',0,...
                               'plotit',true);
    clear a0 c0 h gammaexp c d2 x* y* z numobs n iid glacier a b 