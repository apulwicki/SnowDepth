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

%% Elevation and SD (McGrath style plots of depth)

% Depth vs Elevation (GPS elev)
clf
glacier = 'G02';

glacier_list = ['G04';'G02';'G13']; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)

for j = 1:3;
    glacier = glacier_list(j,:);
    
    index = find(SWE(j).utm(:,3)~=0);
    elev = SWE(j).utm(index,3);
    depth = SWE(j).depth(index,1);
    
    binSize = 10;
    bins = [min(elev):binSize:max(elev)]';
    depthBin = zeros(length(bins)-1,1);
    for i = 1:length(bins)-1
        logiEL = elev>bins(i) & elev<=bins(i+1);
        depthBin(i) = std(depth(logiEL))%/mean(depth(logiEL))*100;
    end

    binsPlot.(glacier) = bins(~isnan(depthBin))-binSize/2;
    depthPlot.(glacier) = depthBin(~isnan(depthBin));

    [myfit.(glacier), gof.(glacier)] = fit(binsPlot.(glacier), depthPlot.(glacier),'poly1');
    p.(glacier) = plot(myfit.(glacier), binsPlot.(glacier), depthPlot.(glacier),'.'); hold on
    set(p.(glacier)(1,1), 'Color',options.RGB(j,:),'MarkerSize', 10)
    set(p.(glacier)(2,1), 'Color',options.RGB(j,:), 'LineWidth', 2)
end

%Elevation span lines
plot([nanmin(nanmin(topo_full_ns.G4.elevation)) nanmax(nanmax(topo_full_ns.G4.elevation))], [68 68], 'Color',options.RGB(1,:), 'LineWidth', 4);
plot([nanmin(nanmin(topo_full_ns.G2.elevation)) nanmax(nanmax(topo_full_ns.G2.elevation))], [66 66], 'Color',options.RGB(2,:), 'LineWidth', 4);
plot([nanmin(nanmin(topo_full_ns.G13.elevation)) nanmax(nanmax(topo_full_ns.G13.elevation))], [64 64], 'Color',options.RGB(3,:), 'LineWidth', 4);
ylim([0 70]); %ylim([0 100]); 
xlim([1900 3100])

xlabel('Elevation (m a.s.l.)'); ylabel([{'Standard deviation as percent'}, {'of mean snow depth (%)'}])
legend([p.G04(1,1) p.G02(1,1) p.G13(1,1)], ...
    {['G04 R^2 = ', num2str(round(gof.G04.rsquare,2))],...
    ['G02 R^2 = ', num2str(round(gof.G02.rsquare,2))],...
    ['G13 R^2 = ', num2str(round(gof.G13.rsquare,2))]}, 'Location', 'east');

    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 8 4.5];
filename = 'binned_std';%_percent';
print([options.path1, filename],'-dpng'); print([options.path2, filename],'-dpng')


%% Normalized Elev

elevG = struct('G04',[1573 2854],'G02',[1906 3098],'G13',[1775 3037]);


for i = 1:3;
    glacier = glacier_list(i,:);
    binsPlot.(glacier)(:,2) = (binsPlot.(glacier)(:,1)- elevG.(glacier)(1,1))/...
        (elevG.(glacier)(1,2)- elevG.(glacier)(1,1));
    
    [myfit.(glacier), gof.(glacier)] = fit(binsPlot.(glacier)(:,2), depthPlot.(glacier),'poly1');
    p.(glacier) = plot(myfit.(glacier),binsPlot.(glacier)(:,2), depthPlot.(glacier),'-'); hold on
end
title('Mean snow depth (binned) vs GPS elevation')
xlabel('Normalized elevation'); ylabel('Mean Snow Depth (cm)')
legend([p.G04(1,1) p.G02(1,1) p.G13(1,1)], ...
    {['G04 R^2 = ', num2str(round(gof.G04.rsquare,2))],...
    ['G02 R^2 = ', num2str(round(gof.G02.rsquare,2))],...
    ['G13 R^2 = ', num2str(round(gof.G13.rsquare,2))]});

%% pulldataSWE test

    %pulldataSWE(data, glacier, pattern, book, person)
z = pulldataSWE(SWE, 'all','all','all','all');
    %scatter(z.utm(:,1),z.utm(:,2), 10, z.depth,'filled');

%% Transect portions variograms -> not enough points....

%Transverse
index(1).transverse = [4:36,49:56,82:92,110:119]+.1;
index(2).transverse = [267:275,239:247,344:354,519:528,371:378,450:458,407:413,483:488]+.2;
index(3).transverse = [783:792,745:760,660:677,628:644,588:602]+.3;

[~, I1] = intersect(SWE(1).label, categorical(index(1).transverse));
[~, I2] = intersect(SWE(2).label, categorical(index(2).transverse));
[~, I3] = intersect(SWE(3).label, categorical(index(3).transverse));

data = [SWE(1).swe(I1), SWE(1).utm(I1,1:2); SWE(2).swe(I2), SWE(2).utm(I2,1:2);...
    SWE(3).swe(I3), SWE(3).utm(I3,1:2)];

d = variogramAlex(data, 20, 400);
vario = variofitAlex(d,'test',1);

%% Normality tests

run OPTIONS.m
options.ZZ = 2; %exclude zigzags
run MAIN

for g = 1:3
glacier = char(options.glacier(g));
pat = categories(SWE(1).pattern);   %pat = pat(1:end-1);

%Overall glacier
%     [~,p,stats] = chi2gof(SWE(g).depth, 'Alpha', 0.1);
%     display([char(glacier),' ', num2str(p), ' ', num2str(stats.chi2stat)])  
    %display([char(glacier),' ', num2str(round(nanstd(SWE(g).depth)/ nanmean(SWE(g).depth)*100)),'%'])
    display([char(glacier),' ', num2str(round(nanmean(SWE(g).depth),1))])

    %Pattern and Person
    per     = categories(SWE(1).person); per = per(1:4);    person = [];
    for i = 1:length(pat)
%         [~,p,stats] = chi2gof(SWE(g).depth(SWE(g).pattern==pat(i)), 'Alpha', 0.1);
%         pattern(i) = stats.chi2stat;
%         pattern(i) = round(nanstd(SWE(g).depth(SWE(g).pattern==pat(i)))/...
%             nanmean(SWE(g).depth(SWE(g).pattern==pat(i)))*100);
        pattern(i) = round(nanmean(SWE(g).depth(SWE(g).pattern==pat(i))),1);
       for k = 1:4
            data = SWE(g).pattern==pat(i) & SWE(g).person==per(k);
%         	[~,p,stats] = chi2gof(SWE(g).depth(data), 'Alpha', 0.1);
%             person(i,k) = [stats.chi2stat];
%             person(i,k) = round(nanstd(SWE(g).depth(data))/...
%                 nanmean(SWE(g).depth(data))*100);
            person(i,k) = round(nanmean(SWE(g).depth(data)),2);
        end
    end
pointvar.(glacier) = table(pattern', 'RowNames',pat);
pointvar.(glacier) = [pointvar.(glacier), ...
    table(person(:,1),person(:,2),person(:,3),person(:,4), 'VariableNames',per, 'RowNames',pat)];
end
