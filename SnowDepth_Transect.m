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

expression = {'stream?','rock'};
%expression = {'None'};
SDcomments = commentsearch(z, expression, 'in');
    
%% Std in vs out channel
    %pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all','all','all','all',1,'fat');

    %summary(z(5).comments)
expression = {'Channel area','Channel?','In channel','stream?','Channel','Probably channel','Channel (320+)',...
    '? Not ice!','? no ice','Not ice','Not on ice!'};
SD_outchannel = commentsearch(z, expression, 'out');
SD_inchannel = commentsearch(z, expression, 'in');

display(['Average std in channel = ', num2str(mean(SD_inchannel(:,3)))])
display(['Average std out of channel = ', num2str(mean(SD_outchannel(:,3)))])


figure(1)
errorbar(SD_outchannel(:,1), SD_outchannel(:,2), SD_outchannel(:,3),'o')
    xlabel('Waypoint number')
    ylabel('Mean snowdepth (cm)')    

figure(2)
scatter(SD_outchannel(:,1), SD_outchannel(:,3))
    xlabel('Waypoint number')
    ylabel('Standard deviation')    
%% Variogram - transect

glacier = 'G04'; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
z1 = pulldata(SD,'Extra',glacier,'Extra','Extra',1,'fat'); %ExtraSD data from nontransect measurements

x = [z(5).depth(:,6);z1(5).depth(:,42)]; x2 = nanmax(x)-x; %convert easting to distance in m
y = [z(5).depth(:,7);z1(5).depth(:,43)]; y2 = nanmax(y)-y; %convert easting to distance in m
z = [nanmean(z(5).depth(:,1:4),2);nanmean(z1(5).depth(:,1:40),2)];

% if ishandle(f1) %clears data from open plots
%     clf(f1); clf(f2);
% end
% 
% %Variogram
% f1 = figure(1); 
% subplot(2,2,1)
%     scatter(x2,y2,4,z,'filled'); box on;
%     ylabel('y'); xlabel('x'); 
%     c = colorbar; c.Label.String = 'Snow depth (cm)';
%     title(glacier)
% subplot(2,2,2)
%     hist(z,20)
%     ylabel('frequency'); xlabel('z')
%     title('histogram of z-values')
% subplot(2,2,3)
%     [d, iid]= variogram([x2 y2],z,'plotit',true,'nrbins',100);
%     title('Isotropic variogram')
% subplot(2,2,4)
%     d2 = variogram([x2 y2],z,'plotit',true,'nrbins',100,'anisotropy',true);
%     title('Anisotropic variogram')
% 
% %variogram fit
% f2 = figure(2);
% subplot(2,1,2)
%     hist(iid(:,3),100)
%     ylabel('frequency'); xlabel('binned lag distance');
% subplot(2,1,1)
%     h=d.distance;
%     gammaexp = d.val;
%     numobs = d.num;
%     a0 = 500; % initial value: range 
%     c0 = 500; % initial value: sill     
%     [a,b,n] = variogramfit(h,gammaexp,a0,c0,numobs,...
%                            'solver','fminsearchbnd',...
%                            'nugget',0,...
%                            'plotit',true);
% clear a0 c0 h gammaexp c d2
 %Variogram
    figure 
    subplot(2,1,1)
        scatter(x2,y2,4,z,'filled'); box on;
        ylabel('y (m)'); xlabel('x (m)'); 
        c = colorbar; c.Label.String = 'Snow depth (cm)';
        title(glacier)
    subplot(2,1,2)
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

%Normality test (One-sample Kolmogorov-Smirnov test)
if kstest(z)
    display('Depth is normally distributed');
else
    display('Depth is NOT normally distributed');
end


%% Variogram for different sections (lower, upper)

%Upper
    glacier = 'G04'; %select data from chosen glacier
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

    [d, ~]= variogram([x2 y2],z,'plotit',false,'nrbins',100);

    %Variogram
    figure 
    subplot(2,1,1)
        scatter(x2,y2,4,z,'filled'); box on;
        ylabel('y (m)'); xlabel('x (m)'); 
        c = colorbar; c.Label.String = 'Snow depth (cm)';
        title({glacier,'Upper Ablation'})
    subplot(2,1,2)
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


%% Autocorrelation
    %Can only do for UT on G4 because it is continuous and evenly spaced data
    %Also for just one SD because they were taken every 60 m (for the most
    %part)
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'SD1','G02','all','UM',1,'fat'); %transect data  
SDmean = nanmean(z(5).depth(:,1:4),2);

autocorr(SDmean) %do autocorrelation of the data *change total lags with size of matrix
%autocorr(SDmean,floor(size(SDmean,1)/4)) %do autocorrelation of the data *change total lags with size of matrix


%% Stats between categories

z = pulldata(SD,'all','G04','all','all',1,'skinny'); 

%Mean SD and variance for each group
[xbar,s2,grp] = grpstats(z(2).depth,z(2).pattern,{'mean','var','gname'})

%One-way ANOVA
[~,~,stats] = anova1(z(2).depth,z(2).pattern)
[c,~,~,gnames] = multcompare(stats);
[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))] %diplays: groups compared, lower CI limit, difference between means, upper CI, p

%n-way ANOVA
p = anovan(z(2).depth, {z(2).book, z(2).person, z(2).pattern},'varnames',{'book','person','pattern'})
%n-way ANOVA with two-factor interactions -> DNE
p = anovan(z(2).depth, {z(2).book, z(2).person, z(2).pattern, z(2).glacier}, ...
    'model','interaction','varnames',{'book','person','pattern','glacier'})

%% Std Stats

glacier = 'G13'; %select data from chosen glacier
pattern = 'UT';
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all',glacier,'all',pattern,1,'fat'); %transect data  

%One-way ANOVA
[~,~,stats] = anova1(nanstd(z(5).depth(:,1:4),1,2),z(5).person);
[c,~,~,gnames] = multcompare(stats);
[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))] %diplays: groups compared, lower CI limit, difference between means, upper CI, p


%% FFT on transects

%z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'SD2','G02','all','UM',1,'fat'); %transect data  
SDmean = nanmean(z(5).depth(:,1:4),2);

Y = fft(SDmean);

Fs = 30;            % Sampling frequency
T = 1/Fs;             % Sampling period
L = size(Y,1);             % Length of signal
t = (0:L-1)*T;        % Time vector

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')