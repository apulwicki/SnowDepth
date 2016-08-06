function z = pulldata(data, book, glacier, person, pattern, quality, format)
% Used for getting subsets of transect data based on categorical values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %This function produces a structued array (z) that contains the selected
    %data, which is done by specifiying book, glaicer, person, pattern, and
    %quality when the function is called. The format of the 'filtered' cell 
    %can be specificed as either 'fat' or 'skinny'. 
    
    % When format 'fat' is selected, the structured array has rows 1-4 
    %with SD1, SD2, SD3, and ExtraSD depth cells, in original formating, 
    %with desired values retained and NaN for all other values. The last row
    %has the 'filtered' depth values along with WP and WP coordinates. The 
    %first column has depth values and the second has comments associated
    %with with WP
        %This format is useful for calculating mean, std, etc.
    %When format 'skinny' is selected, the data is arranged in one column
    %with the measurements (D1, D2, D3, D4) stacked and then the books
    %stacked on eachother as well. The first row of z is the full data with
    %NaN where data is not call and the second row of z is the 'filtered'
    %data. The structure z has columns that correspond to the original data
    %but also arranged as nx1 vectors
        %This format is useful for statistical comparison of groups

    %***Note that you cannot get transect and zigzag data together, must run
    %two pulldata function with transect and zigzag values separately
    
    %       Alexandra Pulwicki July 2016
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set categories if 'all' is selected for glacier, person, or pattern
    if strcmp(glacier,'all')
        glacier = {'G13','G02','G04'};
    end

    if strcmp(person,'all')
       person =  {'AP','GF','CA','AC'};
    end

    if strcmp(pattern,'all')
        pattern = {'B','BT','LC','LH','LM','UC','UH','UM','UT'};
    end

% Set which SD books to look through
    %Range corresponds to the row of the data structure to go through
    %All other snowdepth book data are set to NaNs 
    if strcmp(book,'all')
        range = 1:size(data,2)-1;
        display('This does not include "Extra" data')
    elseif strcmp(book,'SD1')
        range = 1;
        data(2).depth = nan(size(data(2).depth));
        data(3).depth = nan(size(data(3).depth));
        data(4).depth = nan(size(data(4).depth));
        display('This does not include "Extra" data')
    elseif strcmp(book,'SD2') 
        range = 2;
        data(1).depth = nan(size(data(1).depth));
        data(3).depth = nan(size(data(3).depth));
        data(4).depth = nan(size(data(4).depth));
        display('This does not include "Extra" data')
    elseif strcmp(book,'SD3')
        range = 3;
        data(2).depth = nan(size(data(2).depth));
        data(1).depth = nan(size(data(1).depth));
        data(4).depth = nan(size(data(4).depth));
        display('This does not include "Extra" data')
    elseif strcmp(book,'Extra')
        range = 4;
        data(2).depth = nan(size(data(2).depth));
        data(1).depth = nan(size(data(1).depth));
        data(3).depth = nan(size(data(3).depth));
    end

% Replace values that are not desired with NaN but keep original formating 
for j = range %look in the desired books (set by range value above)
    for i = 1:size(data(j).depth,1) %goes through whole matrix
         if ~ismember(data(j).glacier(i,1),glacier) || ... %if it is not part of the desired glacier, person, or patern, relace with NaN
                    ~ismember(data(j).person(i,1),person) ||...
                    ~ismember(data(j).pattern(i,1),pattern); 
             data(j).depth(i,:) = nan;
         end
    end
    
    if strcmp(quality,'all') %selecting based on quality
        continue %no change is you want all the quality)
    elseif size(data(j).Q,2) == 1 %skips over the ExtraSD data because it is all Q=1
        continue
    elseif quality==1 %chosing good quality data
        data(j).depth(data(j).Q(:,1)=='0',1)=nan; %replace values with Q=0 with NaN
        data(j).depth(data(j).Q(:,2)=='0',2)=nan;
        data(j).depth(data(j).Q(:,3)=='0',3)=nan;
        data(j).depth(data(j).Q(:,4)=='0',4)=nan;
    elseif quality==0 %chosing bad quality data
        data(j).depth(data(j).Q(:,1)=='1',1)=nan; %replace values with Q=1 with NaN
        data(j).depth(data(j).Q(:,2)=='1',2)=nan;
        data(j).depth(data(j).Q(:,3)=='1',3)=nan;
        data(j).depth(data(j).Q(:,4)=='1',4)=nan;
    end
    
end

%Adding decimal to waypoint to indicate book
for i = 1:3 
    data(i).depth(:,5) = data(i).depth(:,5) + 0.1*i;
end

%Filtering and formating data 
if strcmp(format,'fat') %fat format = four depth columns, WP, and WP coordinates
    if strcmp(book,'Extra') %ExtraSD data will have extra columns (which is why you can't call SD and ExtraSD values together with doing some sort of processing) 
        filtered = [data(4).depth(:,2:42), data(4).depth(:,45:46)]; %filtered data
        filteredcomments = data(4).comments; %filtered comments (needed for search comments)
    elseif strcmp(book,'all')
        filtered = [data(1).depth; data(2).depth; data(3).depth]; %stacks SD1,2,3 into one matrix
        filteredcomments = [data(1).comments; data(2).comments; data(3).comments]; %stacks the comments
    else 
        filtered = data(range).depth; %otherwise, just keep the data row you want (specificed by range)
        filteredcomments = data(range).comments; %and the comments you want too
    end
    filteredcomments(all(isnan(filtered(:,1:4)),2),:) = []; %removes all the comments that correspond to not desired values (seen as rows of NaN in the edited matrix)
    filtered(all(isnan(filtered(:,1:4)),2),:) = []; %removes all the not desired values and creates the final filtered data matrix
    
    %Combine the SD matrices and filtered matrix into the final structure
    %of z
    z = struct('depth',{data(1).depth, data(2).depth, data(3).depth, data(4).depth, filtered},...
                'comments', {data(1).comments, data(2).comments, data(3).comments, data(4).comments, filteredcomments});
    
elseif strcmp(format,'skinny') %skinny format = stacked data, one row not filtered, one row filtered
    %Compile all data into vectors (nx1)
    f1 = 'depth';    v1 =  [data(1).depth(:,1); data(1).depth(:,2); data(1).depth(:,3); data(1).depth(:,4);
                            data(2).depth(:,1); data(2).depth(:,2); data(2).depth(:,3); data(2).depth(:,4);
                            data(3).depth(:,1); data(3).depth(:,2); data(3).depth(:,3); data(3).depth(:,4)];
    f2 = 'WP';       v2 =  [data(1).depth(:,5); data(1).depth(:,5); data(1).depth(:,5); data(1).depth(:,5);
                            data(2).depth(:,5); data(2).depth(:,5); data(2).depth(:,5); data(2).depth(:,5);
                            data(3).depth(:,5); data(3).depth(:,5); data(3).depth(:,5); data(3).depth(:,5)];
    f3 = 'book';     v3 =  [data(1).book(:,1); data(1).book(:,1); data(1).book(:,1); data(1).book(:,1);
                            data(2).book(:,1); data(2).book(:,1); data(2).book(:,1); data(2).book(:,1);
                            data(3).book(:,1); data(3).book(:,1); data(3).book(:,1); data(3).book(:,1)];
    f4 = 'glacier';  v4 =  [data(1).glacier(:,1); data(1).glacier(:,1); data(1).glacier(:,1); data(1).glacier(:,1);
                            data(2).glacier(:,1); data(2).glacier(:,1); data(2).glacier(:,1); data(2).glacier(:,1);
                            data(3).glacier(:,1); data(3).glacier(:,1); data(3).glacier(:,1); data(3).glacier(:,1)];
    f5 = 'person';   v5 =  [data(1).person(:,1); data(1).person(:,1); data(1).person(:,1); data(1).person(:,1);
                            data(2).person(:,1); data(2).person(:,1); data(2).person(:,1); data(2).person(:,1);
                            data(3).person(:,1); data(3).person(:,1); data(3).person(:,1); data(3).person(:,1)];
    f6 = 'pattern';  v6 =  [data(1).pattern(:,1); data(1).pattern(:,1); data(1).pattern(:,1); data(1).pattern(:,1);
                            data(2).pattern(:,1); data(2).pattern(:,1); data(2).pattern(:,1); data(2).pattern(:,1);
                            data(3).pattern(:,1); data(3).pattern(:,1); data(3).pattern(:,1); data(3).pattern(:,1)];
    f7 = 'Q';        v7 =  [data(1).Q(:,1); data(1).Q(:,2); data(1).Q(:,3); data(1).Q(:,4);
                            data(2).Q(:,1); data(2).Q(:,2); data(2).Q(:,3); data(2).Q(:,4);
                            data(3).Q(:,1); data(3).Q(:,2); data(3).Q(:,3); data(3).Q(:,4)];
    f8 = 'comments'; v8 =  [data(1).comments(:,1); data(1).comments(:,1); data(1).comments(:,1); data(1).comments(:,1);
                            data(2).comments(:,1); data(2).comments(:,1); data(2).comments(:,1); data(2).comments(:,1);
                            data(3).comments(:,1); data(3).comments(:,1); data(3).comments(:,1); data(3).comments(:,1)];                    
    %Combine the SD matrices in the first row of z
    z = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6,f7,v7,f8,v8);                    
    
    %Filter the data and add  to z
    filterq = ~isnan(v1(:,1)); %filter the data based on where the depth value (v1) is NaN or not
    filtered = struct(f1,v1(filterq),f2,v2(filterq),f3,v3(filterq),... %create new filtered matrix that is nx1
        f4,v4(filterq),f5,v5(filterq),f6,v6(filterq),f7,v7(filterq),f8,v8(filterq));
    z = [z; filtered]; %combine full and filtered data                
end

end