% Converting all depth values to SWE 
%       This script rearranges the depth matrices to join together all the
%       data for each glacier. The output structure is in order of sampling
%       so SWE(1) = G04, SWE(2) = G04, SWE(3) = G13. It also adds elevation
%       data from the SPOT5 DEM to all measurements (found using QGIS). The
%       last part contains the various options for interpolating density
%       values. The final output calculates SWE based on these density
%       values.
%
%       Inputs:         SD, ZZ, and Density structures
%       Outputs:        SWE structure (SWE)

%       Alexandra Pulwicki  Created: August 2016
%                           Updated: September 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Combine all data together and group by glacier

    glacier_list = ['G04';'G02';'G13']; %for selecting data from chosen glacier
    
for i = 1:length(glacier_list) %go through each glacier
    glacier = glacier_list(i,:);
        %z = pulldata(data, book, glacier, person, pattern, quality, format)
    z = pulldata(SD,'all',glacier,'all','all',1,'fat'); % pull transect data  
    z1 = pulldata(SD,'Extra',glacier,'all','Extra',1,'fat'); %pull ExtraSD data from nontransect measurements

    %Joins data from transects, extra, and zigzag measurements for each
    %glacier
    SWE(i).depth = [nanmean(z(5).depth(:,1:4),2); nanmean(z1(5).depth(:,1:40),2); ... %depth
        cell2mat(ZZ.depth(ZZ.glacier==glacier,5))]; 
        x = [z(5).depth(:,6); z1(5).depth(:,42); cell2mat(ZZ.depth(ZZ.glacier==glacier,3))]; %easting
        y = [z(5).depth(:,7); z1(5).depth(:,43); cell2mat(ZZ.depth(ZZ.glacier==glacier,4))]; %northing
    SWE(i).utm = [x y]; %UTM coordinates
    SWE(i).pattern = [z(5).pattern; z1(5).pattern; ZZ.book(ZZ.glacier==glacier)]; %pattern label
    SWE(i).glacier = [z(5).glacier; z1(5).glacier; ZZ.glacier(ZZ.glacier==glacier)]; %glacier label
        C = cell(size(ZZ.person)); C(:) = {'\_'}; %get underscores for ZZ label
        mergeMtx = [ZZ.text(:,1),C,ZZ.text(:,2),C,ZZ.text(:,3)]; %create ZZ label
        ZZlabel = cellstr(cell2mat(mergeMtx)); 
    SWE(i).label = [categorical(z(5).depth(:,5)); z1(5).person; categorical(ZZlabel(ZZ.glacier==glacier))]; %label label
    SWE(i).person = [z(5).person; z1(5).person; ZZ.person(ZZ.glacier==glacier)]; % persona label
    SWE(i).book = [z(5).book; z1(5).book; ZZ.book(ZZ.glacier==glacier)]; %book label
    SWE(i).comments = [z(5).comments; z1(5).comments; ZZ.comments(ZZ.glacier==glacier)]; %comments
end

        clear C depth glacier* i mergeMtx WP_index z x y z1 ZZlabel remove_index

%% Add DEM elevations of measurements

%Import elevations
    DEMelev = xlsread('TransectMeasurementLocations.xlsx', 'Sheet1','D1:D3935'); %get DEM elevations of all measurements

%Assign elevations to correct points
    SWE(1).utm(:,3) = DEMelev(1:length(SWE(1).utm),1);
        bit = length(SWE(1).utm); %start index
    SWE(2).utm(:,3) = DEMelev(bit+1:bit+length(SWE(2).utm),1);
        bit = length(SWE(1).utm)+length(SWE(2).utm); %start index
    SWE(3).utm(:,3) = DEMelev(bit+1:bit+length(SWE(3).utm),1);

        clear DEMelev bit
%% Density options

if options.DensitySWE == 1 %depth
    %do nothing

elseif options.DensitySWE == 2 %Donjek mean density (uniform)       Snowpit
    density = nanmean(cell2mat(Density.snowpit(:,2))); %mean density from all snowpits  
    for i = 1:length(SWE) %for each glacier
        SWE(i).density(:) = density; %uniform density for all measures
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000; %SWE calculation
    end

elseif options.DensitySWE == 3 %Donjek mean density (uniform)       SWE tube
    density = nanmean(nanmean(cell2mat(Density.tube(:,2:19)))); %mean density from all SWE tubes 
    for i = 1:length(SWE) %for each glacier
        SWE(i).density(:) = density; %uniform density for all measures
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;  %SWE calculation
    end
 
elseif options.DensitySWE == 4 %Glacier mean density (uniform)      Snowpit
    index = [5,7,1,4,8,10]; %start and end indices for snowpit densities 
    count = 1;
    for i = 1:length(SWE) %for each glacier
        density = nanmean(cell2mat(Density.snowpit(index(count):index(count+1),2))); %mean density for each glacier
        SWE(i).density(:) = density; %uniform density for each glacier
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000; %SWE calculation
        count = count + 2;
    end
    
elseif options.DensitySWE == 5 %Glacier mean density (uniform)      SWE tube
    index = [1,7,8,14,15,33]; %start and end indices for swe tube densities
    count = 1;
    for i = 1:length(SWE) %for each glacier
        density = nanmean(nanmean(cell2mat(Density.tube(index(count):index(count+1),2:19))));  %mean density for each glacier
        SWE(i).density(:) = density; %uniform density for each glacier
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000; %SWE calculation
        count = count + 2;
    end
    
elseif options.DensitySWE == 6 %Linear elevation regression (variable)      Snowpit    
    index = [5,7,1,4,8,10]; %start and end indices for snowpit densities 
    count = 1;
    for i = 1:length(SWE) %for each glacier
        y = cell2mat(Density.snowpit(index(count):index(count+1),2)); %density values
        x = cell2mat(Density.snowpit(index(count):index(count+1),5)); %elevation values
        [fitDen, gof] = fit(x,y,'poly1'); %linear fit
        density = fitDen.p1*SWE(i).utm(:,3)+fitDen.p2; %calculate density for all values from fit
        SWE(i).density = zeros(size(SWE(i).depth,1),1); 
        SWE(i).density(:) = density; %density is spatially variable
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000; %SWE calculation
        count = count + 2;
    end    
    
elseif options.DensitySWE == 7 %Linear elevation regression (variable)      SWE tube    
    index = [1,7,8,14,15,31];  %start and end indices for swe tube densities
    count = 1;
    for i = 1:length(SWE) %for each glacier
        y = nanmean(cell2mat(Density.tube(index(count):index(count+1),2:19)),2); %density values
        x = cell2mat(Density.tube(index(count):index(count+1),22)); %elevation values
        [fitDen, gof] = fit(x,y,'poly1'); %linear fit
        density = fitDen.p1*SWE(i).utm(:,3)+fitDen.p2; %calculate density for all values from fit
        SWE(i).density = zeros(size(SWE(i).depth,1),1);
        SWE(i).density(:) = density; %density is spatially variable
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000; %SWE calculation
        count = count + 2;
    end    
    
elseif options.DensitySWE == 8 %Inverse distance weighted mean (variable)   Snowpit    
    for i = 1:length(SWE) %for each glacier
    %distance matrix
        Epit = repmat(cell2mat(Density.snowpit(:,3))',length(SWE(i).utm(:,1)),1); %repeating pit easting
        Npit = repmat(cell2mat(Density.snowpit(:,4))',length(SWE(i).utm(:,1)),1); %repeating pit northing
        Eloc = repmat(SWE(i).utm(:,1),1,length(Density.snowpit(:,4))); %repeating locations easting
        Nloc = repmat(SWE(i).utm(:,2),1,length(Density.snowpit(:,4))); %repeating locations northing
       
        X = Eloc-Epit; Y = Nloc-Npit; %separation between locations and pits
        weight = 1./hypot(X,Y); %weight is 1/distance between location and pit
        
   %density  
        D = repmat(cell2mat(Density.snowpit(:,2))',length(SWE(i).utm(:,1)),1); %repeating density matrix
        density = D.*weight; %calculate weighted density
        density = sum(density,2)./sum(weight,2); %add weighted density and divide by sums of weights (inverse distance weighting)
            SWE(i).density = zeros(size(SWE(i).depth,1),1); %initialize SWE field
        SWE(i).density(:) = density; %assign density values
   %swe
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;  %SWE calculation
    end 
    
elseif options.DensitySWE == 9 %Inverse distance weighted mean (variable)   SWE tube    
    for i = 1:length(SWE)
    %distance matrix
        Epit = repmat(cell2mat(Density.tube(:,20))',length(SWE(i).utm(:,1)),1); %repeating tube easting
        Npit = repmat(cell2mat(Density.tube(:,21))',length(SWE(i).utm(:,1)),1); %repeating tube northing
        Eloc = repmat(SWE(i).utm(:,1),1,length(Density.tube(:,20))); %repeating locations easting
        Nloc = repmat(SWE(i).utm(:,2),1,length(Density.tube(:,20))); %repeating locations northing
       
        X = Eloc-Epit; Y = Nloc-Npit; %separation between locations and tube
        weight = 1./hypot(X,Y); %weight is 1/distance between location and tube
        
   %density  
        D = repmat(nanmean(cell2mat(Density.tube(:,2:19)),2)',length(SWE(i).utm(:,1)),1); %repeating density matrix
        density = D.*weight; %calculate weighted density
        density = sum(density,2)./sum(weight,2); %add weighted density and divide by sums of weights (inverse distance weighting)
            SWE(i).density = zeros(size(SWE(i).depth,1),1); %initialize SWE field
        SWE(i).density(:) = density; %assign density values
   %swe
        SWE(i).swe = SWE(i).depth(:)/100.*SWE(i).density(:)/1000;  %SWE calculation
    end 

end

SWE = orderfields(SWE); %alphabetically order the fields   

clear density count i fitDen gof index j weight x X y Y E* N* D
clear SD ZZ

% plot(SWE(1).density,'.'); hold on
% plot(SWE(2).density,'.'); hold on
% plot(SWE(3).density,'.');
% legend('G04','G02','G13')