%test variogram

    %pulldata(data, book, glacier, person, pattern, quality, format)

glacier = 'G04';
pattern = 'UH';

    z = pulldata(SD,'all', glacier,'all', pattern, 1,'fat');

title = [glacier,' ', pattern];
d = variogramAlex(z, 15, 705, title);



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
    %clear a0 c0 h gammaexp c d2 x* y* z numobs n iid glacier a b 