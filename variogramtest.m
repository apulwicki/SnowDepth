%test variogram

    %pulldata(data, book, glacier, person, pattern, quality, format)

glacier = 'G04';
pattern = 'all';

    z = pulldata(SD,'all', glacier,'all', pattern, 1,'fat');
    x = z(5).depth(:,6); x2 = nanmax(x)-x; %convert easting to distance in m
    y = z(5).depth(:,7); y2 = nanmax(y)-y; %convert easting to distance in m
    z1 = nanmean(z(5).depth(:,1:4),2);
figure(1)
title = [glacier,' ', pattern];
d = variogramAlex([z1 x2 y2], 15, 1100, title);


%% Variogram function

z = pulldata(SD,'all',glacier,'all',pattern,1,'fat');
    x = z(5).depth(:,6); x2 = nanmax(x)-x; %convert easting to distance in m
    y = z(5).depth(:,7); y2 = nanmax(y)-y; %convert easting to distance in m
    z1 = nanmean(z(5).depth(:,1:4),2);

d1 = variogram([x2 y2],z1,'plotit',false,'nrbins',47,'maxdist',705);

figure(2) 
        h=d1.distance;
        gammaexp = d1.val;
        numobs = d1.num;
        a0 = 500; % initial value: range 
        c0 = 500; % initial value: sill     
        [a,c,n,S] = variogramfit(h,gammaexp,a0,c0,[],...
                               'solver','fminsearchbnd',...
                               'nugget',0,'plotit',true,...
                               'model','gaussian');
    clear a* c* h gammaexp d2 x* y* numobs n iid b 
    
    
    %%
    
    x = rand(1000,1)*4-2;  
    y = rand(1000,1)*4-2;
    z = 3*sin(x*15)+ randn(size(x));
    
 d1 = variogram([x y],z,'plotit',false,'nrbins',50);   
 d = variogramAlex([z x y], 0.0565, 2.825, 'Alex');   
    
d1.val(1:5)
d.val(1:5)
    