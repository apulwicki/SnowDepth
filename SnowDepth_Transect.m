%% Working with transect data
%***** add intro
run SnowDepth_Import.m %Imports snow depth and measurement location data

% Getting easting and northing for each data point (columns 6 and 7) from
% the vector 'closest' from the MeasurementLocation.m script
for i = 1:size(SD1,1) %for all the waypoints in the SD1 matrix
    ind = find(SD1(i,5)==floor(closest(:,3))); % get index of corresponding waypoint in 'closest' (there will be three because 4.1, 4.2, 4.3)
    SD1(i,6:7) = closest(ind(1),1:2); %assign the first easting and northing from 'closest' to column 6 and 7 of SD
end

for i = 1:size(SD2,1)
    ind = find(SD2(i,5)==floor(closest(:,3)));
    SD2(i,6:7) = closest(ind(2),1:2);
end

for i = 1:size(SD3,1)
    ind = find(SD3(i,5)==floor(closest(:,3)));
    SD3(i,6:7) = closest(ind(3),1:2);
end

clear ind i 

% Getting all waypoints in SD1, SD2, and SD3 data (empty for ones not 
% measured). This allows SD matrices to be the same size and have aligned
% WPs and indices
for i=1:max(SD1(:,5))-3 %use WP range so that gaps can be filled
   if i+3~=SD1(i,5) %determine if WP is missing from the data (first WP is 4 so each WP# = index+3)
       SD1 = [SD1(1:i-1,:); nan,nan,nan,nan,i+3,nan,nan ;SD1(i:end,:)]; %insert nan row with WP# in SD1
       SD1_Date = [SD1_Date(1:i-1,1); '<undefined>' ;SD1_Date(i:end,1)]; %insert and undefinied into categorical arrays
       SD1_Glacier = [SD1_Glacier(1:i-1,1); '<undefined>' ;SD1_Glacier(i:end,1)];
       SD1_Pattern = [SD1_Pattern(1:i-1,1); '<undefined>' ;SD1_Pattern(i:end,1)];
       SD1_Person = [SD1_Person(1:i-1,1); '<undefined>' ;SD1_Person(i:end,1)];
       SD1_Book = [SD1_Book(1:i-1,1); '<undefined>' ;SD1_Book(i:end,1)];
       SD1_Q = [SD1_Q(1:i-1,1:4);{'<undefined>', '<undefined>','<undefined>','<undefined>'};SD1_Q(i:end,1:4)];
       SD1raw = [SD1raw(1:i-1,:); repmat({NaN},1,15) ;SD1raw(i:end,:)]; %insert nan rows into cell arrays
       SD1text = [SD1text(1:i-1,:); repmat({NaN},1,15) ;SD1text(i:end,:)];
   end
   if i+3~=SD2(i,5)
       SD2 = [SD2(1:i-1,:); nan,nan,nan,nan,i+3,nan,nan ;SD2(i:end,:)];
       SD2_Date = [SD2_Date(1:i-1,1); '<undefined>' ;SD2_Date(i:end,1)];
       SD2_Glacier = [SD2_Glacier(1:i-1,1); '<undefined>' ;SD2_Glacier(i:end,1)];
       SD2_Pattern = [SD2_Pattern(1:i-1,1); '<undefined>' ;SD2_Pattern(i:end,1)];
       SD2_Person = [SD2_Person(1:i-1,1); '<undefined>' ;SD2_Person(i:end,1)];
       SD2_Book = [SD2_Book(1:i-1,1); '<undefined>' ;SD2_Book(i:end,1)];
       SD2_Q = [SD2_Q(1:i-1,1:4);{'<undefined>', '<undefined>','<undefined>','<undefined>'};SD2_Q(i:end,1:4)];
       SD2raw = [SD2raw(1:i-1,:); repmat({NaN},1,15) ;SD2raw(i:end,:)];
       SD2text = [SD2text(1:i-1,:); repmat({NaN},1,15) ;SD2text(i:end,:)];
   end
   if i+3~=SD3(i,5)
       SD3 = [SD3(1:i-1,:); nan,nan,nan,nan,i+3,nan,nan ;SD3(i:end,:)];
       SD3_Date = [SD3_Date(1:i-1,1); '<undefined>' ;SD3_Date(i:end,1)];
       SD3_Glacier = [SD3_Glacier(1:i-1,1); '<undefined>' ;SD3_Glacier(i:end,1)];
       SD3_Pattern = [SD3_Pattern(1:i-1,1); '<undefined>' ;SD3_Pattern(i:end,1)];
       SD3_Person = [SD3_Person(1:i-1,1); '<undefined>' ;SD3_Person(i:end,1)];
       SD3_Book = [SD3_Book(1:i-1,1); '<undefined>' ;SD3_Book(i:end,1)];
       SD3_Q = [SD3_Q(1:i-1,1:4);{'<undefined>', '<undefined>','<undefined>','<undefined>'};SD3_Q(i:end,1:4)];
       SD3raw = [SD3raw(1:i-1,:); repmat({NaN},1,15) ;SD3raw(i:end,:)];
       SD3text = [SD3text(1:i-1,:); repmat({NaN},1,15) ;SD3text(i:end,:)];
   end
   
end
clear i 


field1 = 'raw';         value1 = {SD1raw, SD2raw, SD3raw, ExtraSDraw};
field2 = 'depth';       value2 = {SD1, SD2, SD3, ExtraSD};
field3 = 'glacier';     value3 = {SD1_Glacier, SD2_Glacier, SD3_Glacier, ExtraSD_Glacier};
field4 = 'pattern';     value4 = {SD1_Pattern, SD2_Pattern, SD3_Pattern, ExtraSD_Pattern};
field5 = 'person';      value5 = {SD1_Person, SD2_Person, SD3_Person, ExtraSD_Person};
field6 = 'Q';           value6 = {SD1_Q, SD2_Q, SD3_Q, ExtraSD_Q};
field7 = 'book';        value7 = {SD1_Book, SD2_Book, SD3_Book, ExtraSD_Book};

SD = struct(field1,value1,field2,value2,field3,value3,field4,value4,...
    field5,value5,field6,value6,field7,value7);
clear value* field* SD1* SD2* SD3* ExtraSD*

%% Basic Stats

    

%% Function test

    %pulldata(data, book, glacier, person, pattern, quality, format)
%z = pulldata(SD,'ZZ','G04','ZZ','ZZ',1,'fat');
z = pulldata(SD,'SD1','G04','AP','UT',1,'fat');

SDmean = nanmean(z(5).depth(:,1:4),2);
SDstd = nanstd(z(5).depth(:,1:4),1,2);  %std normalized by n

errorbar(z(5).depth(:,5), SDmean, SDstd)
xlabel('Waypoint number')
ylabel('Mean snowdepth (cm)')
%% Variogram - transect

glacier = 'G13'; %select data from chosen glacier
% x = [SD1(SD1_Glacier==glacier,6);   SD2(SD2_Glacier==glacier,6); ...
%         SD3(SD3_Glacier==glacier,6);    ExtraSD(ExtraSD_Glacier==glacier,46)]; %easting
% x2 = nanmax(x)-x; %convert easting to distance in m
% y = [SD1(SD1_Glacier==glacier,7);   SD2(SD2_Glacier==glacier,7);...
%         SD3(SD3_Glacier==glacier,7);    ExtraSD(ExtraSD_Glacier==glacier,47)]; %northing
% y2 = nanmax(y)-y; %convert easting to distance in m
% z = [SD1Mean(SD1_Glacier==glacier,1); SD2Mean(SD2_Glacier==glacier,1);...
%         SD3Mean(SD3_Glacier==glacier,1);    ExtraSDMean(ExtraSD_Glacier==glacier,1)]; %mean snow dpeth

%z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all','G13','all','all','1','fat');    
x = z(5).depth(:,6); x2 = nanmax(x)-x; %convert easting to distance in m
y = z(5).depth(:,7); y2 = nanmax(y)-y; %convert easting to distance in m
z = nanmean(z(5).depth(:,1:4),2);

% if ishandle(f1) %clears data from open plots
%     clf(f1); clf(f2);
% end

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




