global options

%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SAVING PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %SFU 
%options.path1 = '/home/glaciology1/Documents/Data/Plots/'; %Plots folder
%options.path2 = '/home/glaciology1/Documents/MastersDocuments/Methods/'; %Latex

    %Laptop
options.path1 = '/Users/Alexandra/Documents/SFU/Data/Plots/'; %Plots
options.path2 = '/Users/Alexandra/Documents/SFU/MastersDocuments/Methods/'; %Latex

%%%%%%%%%%%%%%%%%%%%%% DETREND FED. SAMPLER DENSITY %%%%%%%%%%%%%%%%%%%%%%%

options.TubeDensity     = 2;
                % 1 original fed. sampler density values
                % 2 depth detrended fed. sampler density values

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZIGZAG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Measurement locations
options.ZigzagLocation  = 2;
                % 1 calculate location based on vertex utm
                % 2 calculate location based on last probe point
% Converting to SWE
options.ZigzagSWE       = 2;
                % 1 depth value (raw)
                % 2 SWE value using mean SWE tube value at zigzag

                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLACIER SWE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Converting to SWE

options.DensitySWE     = 8;
                % 1 Depth (raw)
                
                % 2 Donjek mean density (uniform)       Snowpit
                % 3 Donjek mean density (uniform)       SWE tube
                % 4 Glacier mean density (uniform)      Snowpit
                % 5 Glacier mean density (uniform)      SWE tube
                
                % 6 Linear elevation regression (variable)      Snowpit
                % 7 Linear elevation regression (variable)      SWE tube
                % 8 Inverse distance weighted mean (variable)   Snowpit
                % 9 Inverse distance weighted mean (variable)   SWE tube

                
                
                
                
                