%% Import data
run SnowDepth_Import.m %Imports snow depth and measurement location data

%% Basic Stats

    %pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'SD1','G04','AP','UT',1,'fat');

SDmean = nanmean(z(5).depth(:,1:4),2);
SDstd = nanstd(z(5).depth(:,1:4),1,2);  %std normalized by n

figure
errorbar(z(5).depth(:,5), SDmean, SDstd)
xlabel('Waypoint number')
ylabel('Mean snowdepth (cm)')    

%% Function test

    %pulldata(data, book, glacier, person, pattern, quality, format)
z1 = pulldata(SD,'Extra','all','Extra','Extra',1,'fat');
z = pulldata(SD,'all','G04','all','all',1,'skinny');


SDmean = [nanmean(z(5).depth(:,1:4),2); nanmean(z1(5).depth(:,1:40),2)];
SDstd = [nanstd(z(5).depth(:,1:4),1,2); nanstd(z1(5).depth(:,1:40),1,2)];  %std normalized by n

%% Searching comments
    %pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'SD1','all','all','all','all','fat');
SDmean = nanmean(z(5).depth(:,1:4),2);
SDstd = nanstd(z(5).depth(:,1:4),1,2);
    %To search comments in all books you need to do them one at a time...
expression = 'stream|channel';
present = regexpi(z(5).comments, expression);
waypoint = z(5).depth(~cellfun(@isempty,present),5);
index = ismember(z(5).depth(:,5),waypoint); 

z(5).depth(index,1:4)
%SDmean(index)
%SDstd(index)

% Display comment and WP#
all_comments = [num2cell(waypoint), z(5).comments(index,1)];

    clear expression present waypoint index
    
%% Std in vs out channel
    %pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all','all','all','all','all','fat');

[fullSDmean, fullSDstd, all_comments] = commentsearch(z, 'stream|channel','in');

figure(1)
errorbar(fullSDmean(:,2), fullSDmean(:,1), fullSDstd(:,1),'o')
    xlabel('Waypoint number')
    ylabel('Mean snowdepth (cm)')    

figure(2)
scatter(fullSDstd(:,2),fullSDstd(:,1))
    xlabel('Waypoint number')
    ylabel('Standard deviation')    
%% Variogram - transect

glacier = 'G13'; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all',glacier,'all','all','1','fat'); %transect data  
z1 = pulldata(SD,'Extra',glacier,'Extra','Extra',1,'fat'); %ExtraSD data from nontransect measurements

x = [z(5).depth(:,6);z1(5).depth(:,42)]; x2 = nanmax(x)-x; %convert easting to distance in m
y = [z(5).depth(:,7);z1(5).depth(:,43)]; y2 = nanmax(y)-y; %convert easting to distance in m
z = [nanmean(z(5).depth(:,1:4),2);nanmean(z1(5).depth(:,1:40),2)];

if ishandle(f1) %clears data from open plots
    clf(f1); clf(f2);
end

%Variogram
f1 = figure(1); 
subplot(2,2,1)
    scatter(x2,y2,4,z,'filled'); box on;
    ylabel('y'); xlabel('x'); 
    c = colorbar; c.Label.String = 'Snow depth (cm)';
    title(glacier)
subplot(2,2,2)
    hist(z,20)
    ylabel('frequency'); xlabel('z')
    title('histogram of z-values')
subplot(2,2,3)
    [d, iid]= variogram([x2 y2],z,'plotit',true,'nrbins',100);
    title('Isotropic variogram')
subplot(2,2,4)
    d2 = variogram([x2 y2],z,'plotit',true,'nrbins',100,'anisotropy',true);
    title('Anisotropic variogram')

%variogram fit
f2 = figure(2);
subplot(2,1,2)
    hist(iid(:,3),100)
    ylabel('frequency'); xlabel('binned lag distance');
subplot(2,1,1)
    h=d.distance;
    gammaexp = d.val;
    numobs = d.num;
    a0 = 500; % initial value: range 
    c0 = 500; % initial value: sill     
    [a,b,n] = variogramfit(h,gammaexp,a0,c0,numobs,...
                           'solver','fminsearchbnd',...
                           'nugget',0,...
                           'plotit',true);
clear a0 c0 h gammaexp c d2
 

%Normality test (One-sample Kolmogorov-Smirnov test)
if kstest(z)
    display('Depth is normally distributed');
else
    display('Depth is NOT normally distributed');
end

%% Autocorrelation
    %Can only do for UT on G4 because it is continuous and evenly spaced data
    %Also for just one SD because they were taken every 60 m (for the most
    %part)
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'SD1','G02','all','UM','1','fat'); %transect data  
SDmean = nanmean(z(5).depth(:,1:4),2);

autocorr(SDmean,floor(size(SDmean,1)/4)) %do autocorrelation of the data *change total lags with size of matrix


%% Stats between categories

z1 = pulldata(SD,'all','G02','all','UM','1','skinny'); 
z2 = pulldata(SD,'all','G02','all','LM','1','skinny'); 


%Mean SD and variance for each group
[xbar,s2,grp] = grpstats([z1; z2],group,{'mean','var','gname'})
           
%One-way ANOVA
[~,~,stats] = anova1(allSD,group)
[c,~,~,gnames] = multcompare(stats);
[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))] %diplays: groups compared, lower CI limit, difference between means, upper CI, p




