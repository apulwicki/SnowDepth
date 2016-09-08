% Converting all depth values to SWE 

%% Combine all data together

glacier_list = ['G04';'G02';'G13']; %select data from chosen glacier
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
for i = 1:length(glacier_list)
    
    glacier = glacier_list(i,:);
z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
z1 = pulldata(SD,'Extra',glacier,'all','Extra',1,'fat'); %ExtraSD data from nontransect measurements

SWE(i).depth = [nanmean(z(5).depth(:,1:4),2); nanmean(z1(5).depth(:,1:40),2); cell2mat(ZZ.depth(ZZ.glacier==glacier,5))];
    x = [z(5).depth(:,6); z1(5).depth(:,42); cell2mat(ZZ.depth(ZZ.glacier==glacier,3))]; 
    y = [z(5).depth(:,7); z1(5).depth(:,43); cell2mat(ZZ.depth(ZZ.glacier==glacier,4))];
SWE(i).utm = [x y];
SWE(i).pattern = [z(5).pattern; z1(5).pattern; ZZ.book(ZZ.glacier==glacier)];
SWE(i).glacier = [z(5).glacier; z1(5).glacier; ZZ.glacier(ZZ.glacier==glacier)];
    C = cell(size(ZZ.person));                                 
    C(:) = {'\_'};
    mergeMtx = [ZZ.text(:,1),C,ZZ.text(:,2),C,ZZ.text(:,3)];
    ZZlabel = cellstr(cell2mat(mergeMtx));
SWE(i).label = [categorical(z(5).depth(:,5)); z1(5).person; categorical(ZZlabel(ZZ.glacier==glacier))];

end

clear C depth glacier* i mergeMtx WP_index z x y z1 ZZlabel

%% Add DEM elevations of measurements
load DEMelev

for i = 1:3
    value = str2double(cellstr(SWE(i).label(:,1))); array = DEMelev(:,1);
    [~, match] = ismember(value, array);
    match(match==0) = [];

    array = str2double(cellstr(SWE(i).label(:,1))); value = DEMelev(match,1);
    [~, match2] = ismember(value, array);
    match2(match2==0) = [];

    SWE(i).utm(match2,3) = DEMelev(match,2);

end

clear match* array value DEMelev i
%% Density options

if options.DensitySWE == 1
    %do nothing

elseif options.DensitySWE == 2
    density = nanmean(cell2mat(Density.snowpit(:,2)));  
    for i = 1:length(SWE)
        SWE(i).density = zeros(size(SWE(i).depth,1),1);
        SWE(i).density(:) = density;
        
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
    end
    SWE = orderfields(SWE); 

elseif options.DensitySWE == 3
    density = nanmean(nanmean(cell2mat(Density.tube(:,2:19))));  
    for i = 1:length(SWE)
        SWE(i).density = zeros(size(SWE(i).depth,1),1);
        SWE(i).density(:) = density;
        
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
    end
    SWE = orderfields(SWE);
    
end

clear density