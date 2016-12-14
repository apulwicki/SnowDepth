global options
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SAVING PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %SFU 
options.path1 = '/home/glaciology1/Documents/Data/Plots/'; %Plots folder
options.path2 = '/home/glaciology1/Documents/MastersDocuments/Methods/'; %Latex

    %Laptop
% options.path1 = '/Users/Alexandra/Documents/SFU/Data/Plots/'; %Plots
% options.path2 = '/Users/Alexandra/Documents/SFU/MastersDocuments/Methods/'; %Latex


%%%%%%%%%%%%%%%%%%%%%%%%%%%% TO ZZ OR NOT TO ZZ %%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.ZZ              = 2;
                % 1 all data, including zigzags
                % 2 no zigzags (just transects)
                % 3 only zigzag

%%%%%%%%%%%%%%%%%%%%%% DETREND FED. SAMPLER DENSITY %%%%%%%%%%%%%%%%%%%%%%%

options.TubeDensity     = 1;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLACIER RGB %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Colour scheme for glaciers
options.RGB            = [9, 132, 103; 224, 187, 2; 130, 75, 135]/255;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EASE OF USE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
options.glacier        = {'G4','G2','G13'};
options.topoVarsUnits  = {'{\itd_C} (m)', '{\itz} (m a.s.l)','\alpha (^{\circ})','{\itm}',...
                            '{\itN}','{\kappa_P} (m^{-1})','{\kappa_T} (m^{-1})','{\itSx}'};
options.topoVars       = {'{\itd_C}', '{\itz}','\alpha','{\itm}',...
                            '{\itN}','{\kappa_P}','{\kappa_T}','{\itSx}'};

 
                   
                