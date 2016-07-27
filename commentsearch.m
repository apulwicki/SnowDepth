function [ mean, std, all_comments ] = commentsearch( z, expression, inout )
%Searches through comments
    %This function is used to search through the comments of a chosen data set
    %(z) and then resturn only those data that match the search criteria. z
    %is chosen through the pulldata.m function. The expression is the text
    %that you want to search (use | for OR). inout is either 'in' or 'out'
    %based on whether you want the expression to be included or excluded. 

%Initialize matrices
    mean            = [];
    std             = [];
    all_comments    = [];

%Searching through comments
    for i = 1:size(z,2)-1 %Look through the first four rows of the input data (SD1,2,3 and Extra)
        SDmean = nanmean(z(i).depth(:,1:end-3),2); %calculate the mean of the snow depth values
        SDstd = nanstd(z(i).depth(:,1:end-3),1,2); %calculate the std of the snow depth values

        present = regexpi(z(i).comments, expression); %see if the expression exists in the comments
        
        if strcmp(inout,'in') %Output the depth values that *include* the expression
            waypoint = z(i).depth(~cellfun(@isempty,present),5); %find the waypoints that have the expressions (ie. not empty)
            index = ~cellfun(@isempty,present); %create logical array with whether the expression exists for that row or not 
        elseif strcmp(inout,'out') %Output the depth values that *exclude* the expression
            waypoint = z(i).depth(cellfun(@isempty,present),5);%find the waypoints that DON'T have the expressions (ie. empty)
            index = cellfun(@isempty,present);
        end
        
        mean = [mean; [SDmean(index),z(i).depth(index,5)+i*0.1]]; %appends the calculated mean and waypoint number (decimal is book #)
        std = [std; [SDstd(index),z(i).depth(index,5)+i*0.1]]; %appends the calculated std and waypoint number (decimal is book #)

        % Display comments and WP#
        all_comments = [all_comments; [num2cell(waypoint+i*0.1), z(i).comments(index,1)]]; %appends the waypoint number (decimal is book #) and comment (useful for checking)
    end
end

