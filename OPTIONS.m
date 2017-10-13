global options
%%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOT SAVING PATH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %SFU 
options.path1       = '/home/glaciology1/Documents/Data/Plots/'; %Plots folder
options.path2       = '/home/glaciology1/Documents/MastersDocuments/Thesis/'; %Latex
options.path3       = '/home/glaciology1/Documents/MastersDocuments/Methods/'; %Latex


    %Laptop
% options.path1     = '/Users/Alexandra/Documents/SFU/Data/Plots/'; %Plots
% options.path2     = '/Users/Alexandra/Documents/SFU/MastersDocuments/Thesis/'; %Latex
% options.path3     = '/Users/Alexandra/Documents/SFU/MastersDocuments/Methods/'; %Latex
%  

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

options.DensitySWE     = 2;
                % 1 Depth (raw)
                
                % 2 Donjek mean density (uniform)       Snowpit
                % 3 Donjek mean density (uniform)       SWE tube
                % 4 Glacier mean density (uniform)      Snowpit
                % 5 Glacier mean density (uniform)      SWE tube
                
                % 6 Linear elevation regression (variable)      Snowpit
                % 7 Linear elevation regression (variable)      SWE tube
                % 8 Inverse distance weighted mean (variable)   Snowpit
                % 9 Inverse distance weighted mean (variable)   SWE tube

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OBS. PER CELL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Number of observations in DEM cell

options.ObsPerCell     = 2;
                % 1 all transect measurement
                % 2 average of all observation in a DEM cell
                
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GLACIER MAPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load('Basic.mat')

% Size of glacier maps
options.mapsize        = mapsize;
                      
% Where map is NaN
options.mapNaN         = mapNaN;

% Sampled grid coordinates
options.ENgrid         = ENgrid;

% RIG boundary
options.rig            = rig; 

% E and N for imagesc
options.E = E;
options.N = N; 

clear rig E N ENgrid mapNaN mapsize
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% EASE OF USE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
options.glacier        = {'G4','G2','G13'};
options.GName          = {'Glacier 4','Glacier 2','Glacier 13'};
options.DenOpt         = {'S1','F1','S2','F2','S3','F3','S4','F4'};
options.topoVarsUnits  = { '{\itz} (m a.s.l)','{\itd_C} (m)','\alpha','{\itm} (^{\circ})',...
                            '{\itN}','{\kappa} (m^{-1})','Sx'};
options.topoVars       = { '{\itz}','{\itd_C}','\alpha','{\itm}',...
                            '{\itN}','{\kappa}','Sx'};
options.densityName    = {'S1','F1','S2','F2','S3','F3','S4','F4'};                        

% Colour scheme for glaciers
options.RGB            = [9, 132, 103; ...
                          224, 187, 2; ...
                          130, 75, 135]/255; 
                   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ZZ STD %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

options.zzstd           = [0.027;   0.035;  0.040];
                