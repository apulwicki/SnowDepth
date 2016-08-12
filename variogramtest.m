%test variogram

    %pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all','G04','all','LH',1,'fat');

data = [nanmean(z(5).depth(:,1:4),2), z(5).depth(:,6:7)];

%% Calculate distance between points (lag)

variotest = zeros(size(data,1)*size(data,1),1);
for i = 1:size(data,1)
    z = data(i,1); x = data(i,2); y= data(i,3);
    for j = 1:size(data,1)
        variotest(i*j,1) = round(EuclideanDistance(x,y,data(j,2),data(j,3)),-1);
        variotest(i*j,2) = (z-data(j,1))^2;
    end
end


variodata = zeros(max(variotest(:,1)),2);
for n = 1:max(variotest(:,1))
    index = find(variotest(:,1)==n);
    variodata(n,1:2) = [n, sum(variotest(index,2))/(2*length(index))];
end

scatter(variodata(:,1),variodata(:,2))

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