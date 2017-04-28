%% Variogram - fake data
nugget = 1;
range = 200;
sill = 10;
lags = 1:10:500;

h1 = lags <= range;
h2 = lags > range;

% Vdata.val = nugget + ( sill*( 1.5*(lags/range) - 0.5*(lags/range).^3)); 
% Vdata.val = [Vdata.val(h1==1), (sill+nugget)*h2(h2==1)]';
Vdata.val = nugget + sill*(1-exp(-lags/range));
Vdata.val = Vdata.val' + rand(size(Vdata.val))'*2;

Vdata.meanDist  = lags';
Vdata.binCentre = lags';
Vdata.num       = ones(size(lags))';

figure(1); clf
    Vfit = variofitAlex(Vdata,'Fake Data',2);

%% Variogram - transect

figure(1); clf
for g = 1%:3    
    glacier = options.glacier{g};

variogram.(glacier) = variogramAlex([SWE(g).swe SWE(g).utm(:,1:2)], 15, 'default');
%subplot(1,3,g);
    fit.(glacier) = variofitAlex(variogram.(glacier),glacier,2);
    %saveFIG(['/home/glaciology1/Documents/Data/Plots/variofull',glacier])
end
%% Variogram for different sections (lower, upper)

glacier_list = ['G04';'G02';'G13']; %select data from chosen glacier
%glacier_list = 'G13';
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
for i = 1:length(glacier_list)
        glacier = glacier_list(i,:);

%UPPER
        %z = pulldata(data, book, glacier, person, pattern, quality, format)
    z1 = pulldata(SD,'all',glacier,'all','UH',1,'fat'); %transect data  
    z2 = pulldata(SD,'all',glacier,'all','UC',1,'fat'); %transect data  
    z3 = pulldata(SD,'all',glacier,'all','UM',1,'fat'); %transect data  
    z4 = pulldata(SD,'all',glacier,'all','UT',1,'fat'); %transect data  
    z5 = pulldata(SD,'all',glacier,'all','BT',1,'fat'); %transect data  

    x = [z1(5).depth(:,6);z2(5).depth(:,6);z3(5).depth(:,6);z4(5).depth(:,6);z5(5).depth(:,6)]; 
        x2 = nanmax(x)-x; %convert easting to distance in m
    y = [z1(5).depth(:,7);z2(5).depth(:,7);z3(5).depth(:,7);z4(5).depth(:,7);z5(5).depth(:,7)]; 
        y2 = nanmax(y)-y; %convert easting to distance in m
    z = [nanmean(z1(5).depth(:,1:4),2);nanmean(z2(5).depth(:,1:4),2);nanmean(z3(5).depth(:,1:4),2);...
        nanmean(z4(5).depth(:,1:4),2);nanmean(z5(5).depth(:,1:4),2)];
    
    variogram = variogramAlex([z x2 y2], 15, 'default');
    fit = variofitAlex(variogram,['Upper ', glacier]);

    filename = strcat('/home/glaciology1/Documents/Data/Plots/vario_upper',glacier);
    %filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/vario_upper',glacier);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 9];
    print(filename,'-dpng','-r0')
    clf
    
%LOWER
        %z = pulldata(data, book, glacier, person, pattern, quality, format)
    z1 = pulldata(SD,'all',glacier,'all','LH',1,'fat'); %transect data  
    z2 = pulldata(SD,'all',glacier,'all','LC',1,'fat'); %transect data  
    z3 = pulldata(SD,'all',glacier,'all','LM',1,'fat'); %transect data  

    x = [z1(5).depth(:,6);z2(5).depth(:,6);z3(5).depth(:,6)]; 
        x2 = nanmax(x)-x; %convert easting to distance in m
    y = [z1(5).depth(:,7);z2(5).depth(:,7);z3(5).depth(:,7)]; 
        y2 = nanmax(y)-y; %convert easting to distance in m
    z = [nanmean(z1(5).depth(:,1:4),2);nanmean(z2(5).depth(:,1:4),2);nanmean(z3(5).depth(:,1:4),2)];
    
    variogram = variogramAlex([z x2 y2], 15, 'default');
    fit = variofitAlex(variogram, ['Lower ', glacier]);
    
    filename = strcat('/home/glaciology1/Documents/Data/Plots/vario_lower',glacier);
   % filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/vario_lower',glacier);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 9];
    print(filename,'-dpng','-r0')
    clf
end
%% Variogram - zigzag

GZZlabel = ['G04 Z3A'; 'G04 Z2A'; 'G04 Z5B'; 'G02 Z5C'; 'G02 Z7A';'G02 Z3B'; 'G13 Z7C';'G13 Z4C'; 'G13 Z3B'; 'G13 Z5A'];


for i = 1:size(GZZlabel,1)
    x = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,3));
    y = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,4));
    z = cell2mat(ZZ.depth(ZZ.index(i):ZZ.index(i+1)-1,5));
    x2 = nanmax(x)-x;
    y2 = nanmax(y)-y;
    
    lag = 1.5; maxlag = 40;
    variogram = variogramAlex([z x2 y2], lag, maxlag);
    fit = variofitAlex(variogram, [GZZlabel(i,:) ' (lag=' num2str(lag) ', maxlag=' num2str(maxlag) ')']);      
   
    %filename = strcat('/home/glaciology1/Documents/Data/Plots/variogram',GZZlabel(i,:));
    filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/Zigzag/variogram',GZZlabel(i,:));
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 9];
    print(filename,'-dpng','-r0')
    clf

end
clear x y x2 y2 z i k c pointsize filename GZZlabel

% figure(1)
%     subplot(2,2,1)
%     scatter(x2,y2,4,z,'filled'); box on;
%     ylabel('y'); xlabel('x')
%     %title(GZZ_lab(i,:))
%     subplot(2,2,2)
%     hist(z,20)
%     ylabel('frequency'); xlabel('z')
%     title('histogram of z-values')
%     subplot(2,2,3)
%     d = variogram([x2 y2],z,'plotit',true,'nrbins',100);
%     title('Isotropic variogram')
%     subplot(2,2,4)
%     d2 = variogram([x2 y2],z,'plotit',true,'nrbins',100,'anisotropy',true);
%     title('Anisotropic variogram')           