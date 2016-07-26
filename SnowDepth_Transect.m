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
z = pulldata(SD,'all','G04','all','all',1,'fat');

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
glacier = 'G13';
pattern = 'LM';
z=[];
for i = 1:length(SD1Mean)
    if ismember(SD1_Glacier(i),glacier) && ismember(SD1_Pattern(i),pattern) %find values in desired glacier and pattern
        z = [z; SD1Mean(i)]; %append data into a data matrix for plotting
    end
%     if ismember(SD2_Glacier(i),glacier) && ismember(SD2_Pattern(i),pattern)
%         z = [z; SD2Mean(i)];
%     end
%     if ismember(SD3_Glacier(i),glacier) && ismember(SD3_Pattern(i),pattern)
%         z = [z; SD3Mean(i)];
%    end
end
autocorr(z,size(z,1)/4) %do autocorrelation of the data *change total lags with size of matrix


%% Stats between categories

% %Compile all data into vectors (nx1)
% allSD = [SD1(:,1); SD1(:,2); SD1(:,3); SD1(:,4);...
%             SD2(:,1); SD2(:,2); SD2(:,3); SD2(:,4);...
%             SD3(:,1); SD3(:,2); SD3(:,3); SD3(:,4)];
% allSD_Person = [SD1_Person;SD1_Person;SD1_Person;SD1_Person;...
%                 SD2_Person;SD2_Person;SD2_Person;SD2_Person;...
%                 SD3_Person;SD3_Person;SD3_Person;SD3_Person];
% allSD_Glacier = [SD1_Glacier;SD1_Glacier;SD1_Glacier;SD1_Glacier;...
%                 SD2_Glacier;SD2_Glacier;SD2_Glacier;SD2_Glacier;...
%                 SD3_Glacier;SD3_Glacier;SD3_Glacier;SD3_Glacier];        
% allSD_Pattern = [SD1_Pattern;SD1_Pattern;SD1_Pattern;SD1_Pattern;...
%                 SD2_Pattern;SD2_Pattern;SD2_Pattern;SD2_Pattern;...
%                 SD3_Pattern;SD3_Pattern;SD3_Pattern;SD3_Pattern];    
% allSD_Q = [SD1_Q(:,1);SD1_Q(:,2);SD1_Q(:,3);SD1_Q(:,4);...
%                 SD2_Q(:,1);SD2_Q(:,2);SD2_Q(:,3);SD2_Q(:,4);...
%                 SD3_Q(:,1);SD3_Q(:,2);SD3_Q(:,3);SD3_Q(:,4)];               

group = allSD_Person; %chose category group

%Mean SD and variance for each group
[xbar,s2,grp] = grpstats(allSD,group,{'mean','var','gname'})
           
%One-way ANOVA
[~,~,stats] = anova1(allSD,group)
[c,~,~,gnames] = multcompare(stats);
[gnames(c(:,1)), gnames(c(:,2)), num2cell(c(:,3:6))] %diplays: groups compared, lower CI limit, difference between means, upper CI, p




