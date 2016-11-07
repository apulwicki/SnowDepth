

    
%Remove zigzag values
fields = fieldnames(SWE(1));
for i = 1:3
    if      options.ZZ == 1
        return
    elseif  options.ZZ == 2 %no zigzags
        data = SWE(i).pattern~='ZZ';
    elseif  options.ZZ == 3 %only zigzags
        data = SWE(i).pattern=='ZZ';
    end

    for j = 1:length(fields)
        SWE(i).(char(fields(j))) = SWE(i).(char(fields(j)))(data,:);
        if  options.ZZ == 3 %no zigzags
            SWE(i).ZZ = char(SWE(i).label); 
            SWE(i).ZZ = categorical(cellstr(SWE(i).ZZ(:,1:end-6)));
        end
    end
    SWE(i).label = removecats(SWE(i).label);
end
    clear data i j fields