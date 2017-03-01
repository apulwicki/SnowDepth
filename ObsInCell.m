function [ SWE ] = ObsInCell( SWE )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

same_cell = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/same_cell.csv', 1, 2);

if length(SWE) == 3;
    div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
            length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];

    for g = 1:3
       sameG = same_cell(div(g,1):div(g,2));

    fields = fieldnames(SWE(g));
        for f = 1:length(fields)
                [A1, I]  = sort(sameG);
        param = char(fields(f));
        %sort everyone
        data     = SWE(g).(param)(I,:);

        T = diff(A1)==0;    T1 = [0;T(1:end-1)]; T = any([T,T1],2);
        sameG_not = unique(A1(T));
        for i = length(sameG_not):-1:1
           ind = find(sameG_not(i)==A1);
           if isnumeric(data)
           data(ind(1,1),1)      = mean(data(ind));
                if strcmp(param,'swe')
                    data(ind(1,1),2) = std(data(ind)); 
                    data(data(:,2)==0,2) = NaN;
                end
           end
           data(ind(2:end,1),:)    = [];
           A1(ind(2:end,1))     = [];
        end
        SWE(g).(param) = data;
        end    
    end
    
    for g = 1:3
       SWE(g).cellstd = SWE(g).swe(:,2);
       SWE(g).swe(:,2) = [];
    end

% elseif length(SWE) == 1;
%     temputm = 
%     [~ 
    
end
    
    
end

