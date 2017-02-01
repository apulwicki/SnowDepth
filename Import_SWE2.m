

%% Average observations in one cell

if options.ObsPerCell==2
    
    same_cell = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/same_cell.csv', 1, 2);

    div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
            length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];

    for g = 1:3
       glacier  = char(options.glacier(g));
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
       glacier  = char(options.glacier(g));    
       SWE(g).cellstd = SWE(g).swe(:,2);
       SWE(g).swe(:,2) = [];
    end
end

        clear A1 data div f fields glacier g i I ind param T T1 same*

%% Generates all swe otpions
% 
% run OPTIONS
% 
% for t = 2:9
% run OPTIONS.m
% options.DensitySWE  = t;
% run MAIN
% 
%   for i = 1:3
%     glacier = char(options.glacier(i)); 
%     sweOPT(t).(glacier) = [SWE(i).swe, SWE(i).utm(:,1:2)];
%   end
% end
% 
%     clear i t glacier
%     
% run OPTIONS
% run MAIN
    
    