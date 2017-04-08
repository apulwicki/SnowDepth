
%% Average observations in one cell

run OPTIONS
if options.ObsPerCell==2   
    [ SWE, topo_sampled ] = ObsInCell(SWE, topo_sampled);
end

%% Standardizing variables

%Keep a copy of non-standardized variables for range of plots of params sampled
topo_sampled_ns = topo_sampled;
topo_full_ns    = topo_full;

%Standardizing variables
params = fieldnames(topo_sampled.G4);
for i = 1:3
name = char(options.glacier(i));
    for t = 1:length(params)
    field = char(params(t));
    
    topo_sampled.(name).(field) = (topo_sampled.(name).(field)-...
        mean(topo_sampled.(name).(field)))/std(topo_sampled.(name).(field));
    topo_full.(name).(field) = (topo_full.(name).(field)-...
        mean(topo_sampled_ns.(name).(field)))/std(topo_sampled_ns.(name).(field));
    end
end

    clear i name topo field j param* div glacier G t X1 Y1 min_dist fields distance centreline corner r
    
%% Sort all topo param structures

order = {'elevation','centreD','aspect','slope','northness','curvature','Sx'};
for i = 1:3
   name     = char(options.glacier(i));
   topo_full.(name)         = orderfields(topo_full.(name),order);
   topo_full_ns.(name)      = orderfields(topo_full_ns.(name),order);   
   topo_sampled.(name)      = orderfields(topo_sampled.(name),order);
   topo_sampled_ns.(name)   = orderfields(topo_sampled_ns.(name),order);
end
    clear i name order
    
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
    
    