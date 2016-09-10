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
SWE(i).person = [z(5).person; z1(5).person; ZZ.person(ZZ.glacier==glacier)];
SWE(i).book = [z(5).book; z1(5).book; ZZ.book(ZZ.glacier==glacier)];
SWE(i).comments = [z(5).comments; z1(5).comments; ZZ.comments(ZZ.glacier==glacier)];

end

clear C depth glacier* i mergeMtx WP_index z x y z1 ZZlabel remove_index

%% Add DEM elevations of measurements

DEMelev = xlsread('TransectMeasurementLocations.xlsx', 'Sheet1','D1:D3935');

SWE(1).utm(:,3) = DEMelev(1:length(SWE(1).utm),1);
    bit = length(SWE(1).utm);
SWE(2).utm(:,3) = DEMelev(bit+1:bit+length(SWE(2).utm),1);
    bit = length(SWE(1).utm)+length(SWE(2).utm);
SWE(3).utm(:,3) = DEMelev(bit+1:bit+length(SWE(3).utm),1);

    clear DEMelev bit
%% Density options

if options.DensitySWE == 1 %depth
    %do nothing

elseif options.DensitySWE == 2 %Donjek mean density (uniform)       Snowpit
    density = nanmean(cell2mat(Density.snowpit(:,2)));  
    for i = 1:length(SWE)
        SWE(i).density(:) = density;
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
    end

elseif options.DensitySWE == 3 %Donjek mean density (uniform)       SWE tube
    density = nanmean(nanmean(cell2mat(Density.tube(:,2:19))));  
    for i = 1:length(SWE)
        SWE(i).density(:) = density;
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
    end
 
elseif options.DensitySWE == 4 %Glacier mean density (uniform)      Snowpit
    index = [5,7,1,4,8,10];
    count = 1;
    for i = 1:length(SWE)
        density = nanmean(cell2mat(Density.snowpit(index(count):index(count+1),2)));
        SWE(i).density(:) = density;
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
        count = count + 2;
    end
    
elseif options.DensitySWE == 5 %Glacier mean density (uniform)      SWE tube
    index = [1,7,8,14,15,33];
    count = 1;
    for i = 1:length(SWE)
        density = nanmean(nanmean(cell2mat(Density.tube(index(count):index(count+1),2:19)))); 
        SWE(i).density(:) = density;
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
        count = count + 2;
    end
    
elseif options.DensitySWE == 6 %Linear elevation regression (variable)      Snowpit    
    index = [5,7,1,4,8,10];
    count = 1;
    for i = 1:length(SWE)
        y = cell2mat(Density.snowpit(index(count):index(count+1),2));
        x = cell2mat(Density.snowpit(index(count):index(count+1),5));
        [fitDen, gof] = fit(x,y,'poly1');
        density = fitDen.p1*SWE(i).utm(:,3)+fitDen.p2;
        SWE(i).density = zeros(size(SWE(i).depth,1),1);
        SWE(i).density(:) = density;
        
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
        count = count + 2;
    end    
    
elseif options.DensitySWE == 7 %Linear elevation regression (variable)      SWE tube    
    index = [1,7,8,14,15,31];
    count = 1;
    for i = 1:length(SWE)
        y = nanmean(cell2mat(Density.tube(index(count):index(count+1),2:19)),2);
        x = cell2mat(Density.tube(index(count):index(count+1),22));
        [fitDen, gof] = fit(x,y,'poly1');
        density = fitDen.p1*SWE(i).utm(:,3)+fitDen.p2;
        SWE(i).density = zeros(size(SWE(i).depth,1),1);
        SWE(i).density(:) = density;
        
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
        count = count + 2;
    end    
    
elseif options.DensitySWE == 8 %Inverse distance weighted mean (variable)   Snowpit    
    for i = 1:length(SWE)
    %distance matrix
            X = zeros(length(SWE(i).utm(:,1)),length(Density.snowpit)); Y = X;
        for j = 1:length(Density.snowpit)
            X(:,j) = SWE(i).utm(:,1)-cell2mat(Density.snowpit(j,3));
            Y(:,j) = SWE(i).utm(:,2)-cell2mat(Density.snowpit(j,4));
        end 
        weight = 1./hypot(X,Y);
        
        %density  
        density = zeros(length(SWE(i).utm(:,1)),length(Density.snowpit));
        for j = 1:length(Density.snowpit)
            density(:,j) = cell2mat(Density.snowpit(j,2))*weight(:,j);
        end     
        density = sum(density,2)./sum(weight,2);

        SWE(i).density = zeros(size(SWE(i).depth,1),1);
        SWE(i).density(:) = density;
        
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
    end 
    
elseif options.DensitySWE == 9 %Inverse distance weighted mean (variable)   SWE tube    
    for i = 1:length(SWE)
    %distance matrix
            X = zeros(length(SWE(i).utm(:,1)),length(Density.tube)); Y = X;
        for j = 1:length(Density.tube)
            X(:,j) = SWE(i).utm(:,1)-cell2mat(Density.tube(j,20));
            Y(:,j) = SWE(i).utm(:,2)-cell2mat(Density.tube(j,21));
        end 
        weight = 1./hypot(X,Y);
        
        %density  
        density = zeros(length(SWE(i).utm(:,1)),length(Density.tube));
        for j = 1:length(Density.tube)
            density(:,j) = nanmean(cell2mat(Density.tube(j,2:19)))*weight(:,j);
        end     
        density = sum(density,2)./sum(weight,2);

        SWE(i).density = zeros(size(SWE(i).depth,1),1);
        SWE(i).density(:) = density;
        
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;
    end 

end

SWE = orderfields(SWE);    

clear density count i fitDen gof index j weight x X y Y 
clear SD ZZ

% plot(SWE(1).density,'.'); hold on
% plot(SWE(2).density,'.'); hold on
% plot(SWE(3).density,'.');
% legend('G04','G02','G13')