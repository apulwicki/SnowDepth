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

expression = {'Channel area','Channel?','In channel','stream?','Channel','Probably channel','Channel (320+)',...
    '? Not ice!','? no ice','Not ice','Not on ice!','Crevasse?'};
inout = 'out';

glacier = 'G13';
z = pulldata(SD,'all',glacier,'all','all',1,'fat');
c13 = commentsearch(z, expression, inout); c13 = c13(:,2);
g13 = cell(length(c13),1); g13(:) = {'G13'};

glacier = 'G02';
z = pulldata(SD,'all',glacier,'all','all',1,'fat');
c02 = commentsearch(z, expression, inout); c02 = c02(:,2);
g02 = cell(length(c02),1); g02(:) = {'G02'};

glacier = 'G04';
z = pulldata(SD,'all',glacier,'all','all',1,'fat');
c04 = commentsearch(z, expression, inout); c04 = c04(:,2);
g04 = cell(length(c04),1); g04(:) = {'G04'};

% display(['Average std in channel = ', num2str(mean(SD_inchannel(:,3)))])
% display(['Average std out of channel = ', num2str(mean(SD_outchannel(:,3)))])

boxplot([c04;c02;c13], [g04;g02;g13])
    ylabel('Snow depth (cm)')
    title({'Snow depth variability between glaciers','Include "channel"'})
    
    clear c* g*
% figure(1)
% errorbar(SD_outchannel(:,1), SD_outchannel(:,2), SD_outchannel(:,3),'o')
%     xlabel('Waypoint number')
%     ylabel('Mean snowdepth (cm)')    
% 
% figure(2)
% scatter(SD_outchannel(:,1), SD_outchannel(:,3))
%     xlabel('Waypoint number')
%     ylabel('Standard deviation')   

%% Box plots

expression = {'Channel area','Channel?','In channel','stream?','Channel','Probably channel','Channel (320+)',...
    '? Not ice!','? no ice','Not ice','Not on ice!','Crevasse?'};

glacier = 'G02'; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
    z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
    z1 = pulldata(SD,'Extra',glacier,'Extra','Extra',1,'fat'); %ExtraSD data from nontransect measurements
z02 = [nanmean(z(5).depth(:,1:4),2);nanmean(z1(5).depth(:,1:40),2)];
g02 = [z(5).glacier; z1(5).glacier];

glacier = 'G04'; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
    z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
    z1 = pulldata(SD,'Extra',glacier,'Extra','Extra',1,'fat'); %ExtraSD data from nontransect measurements
z04 = [nanmean(z(5).depth(:,1:4),2);nanmean(z1(5).depth(:,1:40),2)];
g04 = [z(5).glacier; z1(5).glacier];

glacier = 'G13'; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
    z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
    z1 = pulldata(SD,'Extra',glacier,'Extra','Extra',1,'fat'); %ExtraSD data from nontransect measurements
z13 = [nanmean(z(5).depth(:,1:4),2);nanmean(z1(5).depth(:,1:40),2)];
g13 = [z(5).glacier; z1(5).glacier];

boxplot([z04; z02; z13], [g04; g02; g13], 'GroupOrder',{'G04','G02','G13'})
    ylabel('Snow depth (cm)')
    title('Snow depth variability between glaciers')

    %clear z* g* glacier

%% Variogram - transect

glacier = 'G02'; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
z1 = pulldata(SD,'Extra',glacier,'Extra','Extra',1,'fat'); %ExtraSD data from nontransect measurements

x = [z(5).depth(:,6);z1(5).depth(:,42)]; x2 = nanmax(x)-x; %convert easting to distance in m
y = [z(5).depth(:,7);z1(5).depth(:,43)]; y2 = nanmax(y)-y; %convert easting to distance in m
z = [nanmean(z(5).depth(:,1:4),2);nanmean(z1(5).depth(:,1:40),2)];

[d, SS, SG] = variogramAlex([z x2 y2], 15, 'default', glacier);

    display('spherical')
    SS.nugget
    SS.range
    
    display('gaussian')
    SG.nugget
    SG.range

   %filename = strcat('/home/glaciology1/Documents/Data/Plots/variofull',glacier);
    filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/variofull',glacier);
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 8 9];
    print(filename,'-dpng','-r0')
   
   
    
glacier = 'G04'; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'SD1',glacier,'all','all',1,'fat'); %transect data  

x = z(5).depth(1:2:end-1,6); x2 = nanmax(x)-x; %convert easting to distance in m
y = z(5).depth(1:2:end-1,7); y2 = nanmax(y)-y; %convert easting to distance in m
z = nanmean(z(5).depth(1:2:end-1,1:4),2);

d = variogramAlex([z x2 y2], 15, 'default', glacier);    
    
    
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
%     figure 
%     subplot(2,1,1)
%         scatter(x2,y2,4,z,'filled'); box on;
%         ylabel('y (m)'); xlabel('x (m)'); 
%         c = colorbar; c.Label.String = 'Snow depth (cm)';
%         title(glacier)
%     subplot(2,1,2)
%         h=d.distance;
%         gammaexp = d.val;
%         numobs = d.num;
%         a0 = 500; % initial value: range 
%         c0 = 500; % initial value: sill     
%         [a,b,n] = variogramfit(h,gammaexp,a0,c0,numobs,...
%                                'solver','fminsearchbnd',...
%                                'nugget',0,...
%                                'plotit',true);
%     clear a0 c0 h gammaexp c d2 x* y* z numobs n iid glacier a b 

%Normality test (One-sample Kolmogorov-Smirnov test)
if kstest(z)
    display('Depth is normally distributed');
else
    display('Depth is NOT normally distributed');
end


%% Variogram for different sections (lower, upper)

glacier = 'G13'; %select data from chosen glacier

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
    
    [d, SS, SG] = variogramAlex([z x2 y2], 15, 'default', ['Upper ', glacier]);

    display('spherical')
    SS.nugget
    SS.range
    
    display('gaussian')
    SG.nugget
    SG.range
   %filename = strcat('/home/glaciology1/Documents/Data/Plots/vario_upper',glacier);
%     filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/vario_upper',glacier);
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 8 9];
%     print(filename,'-dpng','-r0')
%     clf
    
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
    
    [d, SS, SG] = variogramAlex([z x2 y2], 15, 'default', ['Lower ', glacier]);

    display('spherical')
    SS.nugget
    SS.range
    
    display('gaussian')
    SG.nugget
    SG.range
    
%    %filename = strcat('/home/glaciology1/Documents/Data/Plots/vario_lower',glacier);
%     filename = strcat('/Users/Alexandra/Documents/SFU/Data/Plots/vario_lower',glacier);
%     fig = gcf;
%     fig.PaperUnits = 'inches';
%     fig.PaperPosition = [0 0 8 9];
%     print(filename,'-dpng','-r0')

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

%% Same shape, check if people different

pattern = 'LH';
glacier = 'G04';

%z = pulldata(data, book, glacier, person, pattern, quality, format)
zAP = pulldata(SD,'all',glacier,'AP',pattern,1,'fat'); %transect data 
zAP_mean = [zAP(5).depth(:,5), nanmean(zAP(5).depth(:,1:4),2), nanstd(zAP(5).depth(:,1:4),1,2)];

zCA = pulldata(SD,'all',glacier,'CA',pattern,1,'fat'); %transect data 
zCA_mean = [zCA(5).depth(:,5), nanmean(zCA(5).depth(:,1:4),2), nanstd(zCA(5).depth(:,1:4),1,2)];

zGF = pulldata(SD,'all',glacier,'GF',pattern,1,'fat'); %transect data 
zGF_mean = [zGF(5).depth(:,5), nanmean(zGF(5).depth(:,1:4),2), nanstd(zGF(5).depth(:,1:4),1,2)];

figure(1)
errorbar(zAP_mean(:,1),zAP_mean(:,2),zAP_mean(:,3)); hold on
errorbar(zCA_mean(:,1),zCA_mean(:,2),zCA_mean(:,3)); hold on
errorbar(zGF_mean(:,1),zGF_mean(:,2),zGF_mean(:,3))
    xlabel('Waypoint')
    ylabel('Depth on G02 UH')
    legend('AP','CA','GF')

%One-way ANOVA
z = pulldata(SD,'all',glacier,'all',pattern,1,'skinny'); 
[~,~,stats] = anova1(z(2).depth,z(2).person)
[c,~,~,gnames] = multcompare(stats);
[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))] %diplays: groups compared, lower CI limit, difference between means, upper CI, p

%% Std Stats

glacier = 'G13'; %select data from chosen glacier
pattern = 'UT';
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all',glacier,'all',pattern,1,'fat'); %transect data  

%One-way ANOVA
[~,~,stats] = anova1(nanstd(z(5).depth(:,1:4),1,2),z(5).person);
[c,~,~,gnames] = multcompare(stats);
[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))] %diplays: groups compared, lower CI limit, difference between means, upper CI, p

%% Comparing std values

glacier = 'G04'; %select data from chosen glacier
pattern = 'all';
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all',glacier,'all',pattern,1,'fat'); %transect data  

std = round(mean(nanstd(z(5).depth(:,1:4),1,2)),1);
display('overall',num2str(std))
cats = categories(z(5).person);
for i= 1:size(cats,1)
    index = z(5).person==cats(i);
    std = round(mean(nanstd(z(5).depth(index,1:4),1,2)),1);
    display(cats(i),num2str(std))
end