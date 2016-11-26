% Including ro exculing zigzag values from the SWE structure 
%       This script filters the SWE values based on the SWE.pattern
%       categorical field to include or exclude zigzag values. Change which
%       option to use in the OPTIONS.m script
%
%       Inputs:         OPTIONS.m
%       Outputs:        SWE structure (SWE)

%       Alexandra Pulwicki  Created: November 2016
%                           Updated: November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %Get field values to make sure all fields end up the same length
fields = fieldnames(SWE(1));
for i = 1:3
     %Create logical matrix that tell you whether you want to keep the
     %value or not
    if      options.ZZ == 1 %all data
        return
    elseif  options.ZZ == 2 %no zigzags
        data = SWE(i).pattern~='ZZ';
    elseif  options.ZZ == 3 %only zigzags
        data = SWE(i).pattern=='ZZ';
    end

     %Go through each field and only keep the data that you want
    for j = 1:length(fields)
        SWE(i).(char(fields(j))) = SWE(i).(char(fields(j)))(data,:);
        if  options.ZZ == 3 %no zigzags
            SWE(i).ZZ = char(SWE(i).label); 
            SWE(i).ZZ = categorical(cellstr(SWE(i).ZZ(:,1:end-6)));
        end
    end
     %Remove unused categorical labels
    SWE(i).label = removecats(SWE(i).label);
end
    clear data i j fields