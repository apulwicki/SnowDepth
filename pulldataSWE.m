function z = pulldataSWE(data, glacier, pattern, book, person)
% Used for getting subsets of transect data based on categorical values
    %This function produces a structued array (z) that contains the selected
    %data, which is done by specifiying book, glaicer, person, pattern, and
    %quality when the function is called. The format of the 'filtered' cell 
    %can be specificed as either 'fat' or 'skinny'. 
    
    % When format 'fat' is selected, the structured array has rows 1-4 
    %with SD1, SD2, SD3, and ExtraSD depth cells, in original formating, 
    %with desired values retained and NaN for all other values. The last row
    %has the 'filtered' depth values along with WP and WP coordinates. The 
    %first column has depth values and the remaining have glacier, person, 
    %pattern, and comments associated with the WP.
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
    
    %       Alexandra Pulwicki  Created: July 2016
    %                           Updated: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set which glacier 
    %Range corresponds to the row of the data structure to go through
    %All other snowdepth book data become blank
    if strcmp(glacier,'all')
        range = 1:size(data,2);
    elseif strcmp(glacier,'G04')
        range = 1;
    elseif strcmp(glacier,'G02') 
        range = 2;
    elseif strcmp(glacier,'G13')
        range = 3;
    end

    index = 1:size(data,2);
    index = ~ismember(index,range);
    name = fieldnames(data);
    for i = find(index)
        for j = 1:length(name)
            data(i).(char(name(j))) = [];
        end
    end
    for j = 1:length(name)
        dataTmp.(char(name(j))) = [data(1).(char(name(j))); data(2).(char(name(j))); data(3).(char(name(j)))];
    end
    data = dataTmp;
    
    % Set categories if 'all' is selected for book, person, or pattern
    if strcmp(book,'all')
        book = categories(data(1).book)';
    end

    if strcmp(person,'all')
       person =  categories(data(1).person)';
    end

    if strcmp(pattern,'all')
        pattern = categories(data(1).pattern)';
    end
    
    index_pattern = ismember(data.pattern,pattern);
    index_book = ismember(data.book,book);
    index_person = ismember(data.pattern,book);
    
    
    
%Filtering and formating data 
    if strcmp(book,'Extra') %ExtraSD data will have extra columns (which is why you can't call SD and ExtraSD values together with doing some sort of processing) 
        filtered = [data(4).depth(:,1:42), data(4).depth(:,44:45)]; %filtered data
        filteredcomments = data(4).comments; %filtered comments (needed for search comments)
        filteredglacier = data(4).glacier; %filtered glacier (needed for search glacier)
        filteredperson = data(4).person; %filtered person (needed for search person)
        filteredpattern = data(4).pattern; %filtered pattern (needed for search pattern)
    elseif strcmp(book,'all')
        filtered = [data(1).depth; data(2).depth; data(3).depth]; %stacks SD1,2,3 into one matrix
        filteredcomments = [data(1).comments; data(2).comments; data(3).comments]; %stacks the comments
        filteredglacier = [data(1).glacier; data(2).glacier; data(3).glacier]; %stacks the glacier
        filteredperson = [data(1).person; data(2).person; data(3).person]; %stacks the person
        filteredpattern = [data(1).pattern; data(2).pattern; data(3).pattern]; %stacks the pattern
    else 
        filtered = data(range).depth; %otherwise, just keep the data row you want (specificed by range)
        filteredcomments = data(range).comments; %and the comments you want too
        filteredglacier = data(range).glacier; %and the glacier you want too
        filteredperson = data(range).person; %and the person you want too
        filteredpattern = data(range).pattern; %and the pattern you want too
    end
    filteredcomments(all(isnan(filtered(:,1:4)),2),:) = []; %removes all the comments that correspond to not desired values (seen as rows of NaN in the edited matrix)
    filteredglacier(all(isnan(filtered(:,1:4)),2),:) = []; %removes all the glacier that correspond to not desired values (seen as rows of NaN in the edited matrix)
    filteredperson(all(isnan(filtered(:,1:4)),2),:) = []; %removes all the person that correspond to not desired values (seen as rows of NaN in the edited matrix)
    filteredpattern(all(isnan(filtered(:,1:4)),2),:) = []; %removes all the pattern that correspond to not desired values (seen as rows of NaN in the edited matrix)
    filtered(all(isnan(filtered(:,1:4)),2),:) = []; %removes all the not desired values and creates the final filtered data matrix
    
    %Combine the SD matrices and filtered matrix into the final structure
    %of z
    f1 = 'depth';    v1 = {data(1).depth, data(2).depth, data(3).depth, data(4).depth, filtered};
    f4 = 'glacier';  v4 = {data(1).glacier, data(2).glacier, data(3).glacier, data(4).glacier, filteredglacier};
    f5 = 'person';   v5 = {data(1).person, data(2).person, data(3).person, data(4).person, filteredperson};
    f6 = 'pattern';  v6 = {data(1).pattern, data(2).pattern, data(3).pattern, data(4).pattern, filteredpattern};
    f8 = 'comments'; v8 = {data(1).comments, data(2).comments, data(3).comments, data(4).comments, filteredcomments};                 
    %Combine the SD matrices in the first row of z
    z = struct(f1,v1,f4,v4,f5,v5,f6,v6,f8,v8);         
            
end