function z = pulldataSWE(data, glacier, pattern, book, person)
% Used for getting subsets of all data based on categorical values
    %This function produces a structued array (z) that contains the selected
    %data, which is done by specifiying book, glaicer, person, and pattern 
    %when the function is called. The resulting structure combines all the
    %data into each field. 
    
    %       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Select glacier
    %Index corresponds to the row of the data structure to *delete*
    %All other glacier data becomes blank
    if strcmp(glacier,'all')
        index = [0 0 0];
    elseif strcmp(glacier,'G04')
        index = [0 1 1];
    elseif strcmp(glacier,'G02') 
        index = [1 0 1];
    elseif strcmp(glacier,'G13')
        index = [1 1 0];
    end
    
    %Deletes the data from other glacier
    name = fieldnames(data); %gets field names to cycle through
    for i = find(index) %find which glacier to delete
        for j = 1:length(name) %cycle through field names
            data(i).(char(name(j))) = []; %make all field names for glaciers to delete blank
        end
    end
    
%Compiles the data from all three glaciers into single columns
    for j = 1:length(name)
        dataTmp.(char(name(j))) = [data(1).(char(name(j))); data(2).(char(name(j))); data(3).(char(name(j)))];
    end
    data = dataTmp;

%Keep only desired data
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
    
    %Get logical index for the data that is desired for each field
    index2 = [ismember(data.pattern,pattern), ismember(data.book,book), ismember(data.person,person)];
    index2 = all(index2,2); %return index for data that falls into call selected field categories
   
    %Select data based on index
    for j = 1:length(name)
        data.(char(name(j))) = data.(char(name(j)))(index2,:);
    end

%Return data matrix
    z = data;
end