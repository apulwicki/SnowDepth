function z = pulldata(data, book, glacier, person, pattern, quality, format)

    %Note that you cannot get transect and zigzag data together, must run
    %two pulldata function with transect and zigzagvalues separately


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
if strcmp(book,'all')
    range = 1:size(data,2)-1;
elseif strcmp(book,'SD1')
    range = 1;
    data(2).depth = nan(size(data(2).depth));
    data(3).depth = nan(size(data(3).depth));
    data(4).depth = nan(size(data(4).depth));
elseif strcmp(book,'SD2') 
    range = 2;
    data(1).depth = nan(size(data(1).depth));
    data(3).depth = nan(size(data(3).depth));
    data(4).depth = nan(size(data(4).depth));
elseif strcmp(book,'SD3')
    range = 3;
    data(2).depth = nan(size(data(2).depth));
    data(1).depth = nan(size(data(1).depth));
    data(4).depth = nan(size(data(4).depth));
elseif strcmp(book,'ZZ')
    range = 4;
    data(2).depth = nan(size(data(2).depth));
    data(1).depth = nan(size(data(1).depth));
    data(3).depth = nan(size(data(3).depth));
end

% Make a new matrix with desired data
for j = range;
    for i = 1:size(data(j).depth,1)
         if ~ismember(data(j).glacier(i,1),glacier) || ...
                    ~ismember(data(j).person(i,1),person) ||...
                    ~ismember(data(j).pattern(i,1),pattern); 
             data(j).depth(i,:) = nan;
         end
    end
    
    if strcmp(quality,'all')
        continue
    elseif size(data(j).Q,2) == 1
        continue
    elseif quality==1
        data(j).depth(data(j).Q(:,1)=='0',1)=nan; 
        data(j).depth(data(j).Q(:,2)=='0',2)=nan;
        data(j).depth(data(j).Q(:,3)=='0',3)=nan;
        data(j).depth(data(j).Q(:,4)=='0',4)=nan;
    elseif quality==0
        data(j).depth(data(j).Q(:,1)=='1',1)=nan; 
        data(j).depth(data(j).Q(:,2)=='1',2)=nan;
        data(j).depth(data(j).Q(:,3)=='1',3)=nan;
        data(j).depth(data(j).Q(:,4)=='1',4)=nan;
    end
    
end

if strcmp(format,'fat')
    if strcmp(book,'ZZ')
            filtered = [data(4).depth(:,2:42), data(4).depth(:,46:47)];
            filtered(all(isnan(filtered),2),:) = [];
    elseif strcmp(book,'all')
        filtered = [data(1).depth; data(2).depth; data(3).depth];
    else 
        filtered = data(range).depth;
    end
    filtered(all(isnan(filtered(:,1:4)),2),:) = [];
    z = struct('depth',{data(1).depth, data(2).depth, data(3).depth, data(4).depth, filtered});
    
elseif strcmp(format,'skinny') %Compile all data into vectors (nx1)
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
    f7 = 'Q';        v7 =  [data(1).Q(:,1); data(1).Q(:,1); data(1).Q(:,1); data(1).Q(:,1);
                            data(2).Q(:,1); data(2).Q(:,1); data(2).Q(:,1); data(2).Q(:,1);
                            data(3).Q(:,1); data(3).Q(:,1); data(3).Q(:,1); data(3).Q(:,1)];
                            
    z = struct(f1,v1,f2,v2,f3,v3,f4,v4,f5,v5,f6,v6,f7,v7);                    
    
    filterq = ~isnan(v1);
    filtered = struct(f1,v1(filterq),f2,v2(filterq),f3,v3(filterq),...
        f4,v4(filterq),f5,v5(filterq),f6,v6(filterq),f7,v7(filterq));
    
    z = [z; filtered];                 
end

end