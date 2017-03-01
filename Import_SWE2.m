

%% Add cell # of observation

    same_cell = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/same_cell.csv', 1, 2);

    div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
            length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];

    for g = 1:3
       SWE(g).cellN = same_cell(div(g,1):div(g,2));
    end

%% Average observations in one cell

if options.ObsPerCell==2   
    SWE = ObsInCell(SWE);
end

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
    
    