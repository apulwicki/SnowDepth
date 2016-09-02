% Converting all depth values to SWE 

%% Combine all data together


glacier_list = ['G04';'G02';'G13']; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
for i = 1:length(glacier_list)
    
    glacier = glacier_list(i,:);
z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
z1 = pulldata(SD,'Extra',glacier,'all','Extra',1,'fat'); %ExtraSD data from nontransect measurements

x = [z(5).depth(:,6); z1(5).depth(:,42); cell2mat(ZZ.depth(ZZ.glacier==glacier,3))]; 
y = [z(5).depth(:,7); z1(5).depth(:,43); cell2mat(ZZ.depth(ZZ.glacier==glacier,4))];
depth = [nanmean(z(5).depth(:,1:4),2); nanmean(z1(5).depth(:,1:40),2); cell2mat(ZZ.depth(ZZ.glacier==glacier,5))];

SWE(i).swe = [depth x y];
SWE(i).pattern = [z(5).pattern; z1(5).pattern; ZZ.book(ZZ.glacier==glacier)];
SWE(i).label = [z(5).depth(:,5); z1(5).person; 
end
