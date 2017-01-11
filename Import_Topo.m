% Importing Topographic Data
%       This script imports topographic data that was exported from QGIS.
%       It imports the DEMs (ratser) as well as the values of topographic
%       parameters at each sampling location (found in QGIS). It also
%       creates a standardized set of params (to be used in MLR for
%       relevant coefficients) as well as a true value set (for plotting).
%
%       Inputs:         GlacierTopos/*, 'SPOT_TopoParams_noZZ.xlsx'
%       Outputs:        Snowdepth structure (SD)
%                       Elevations from GPS WPs (gps_elev)

%       Alexandra Pulwicki  Created: September 2016
%                           Updated: November 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Importing DEMs
run OPTIONS

files   = dir('/home/glaciology1/Documents/Data/GlacierTopos/*.tif'); 
v       = cell(0,0);
Index9999   = [1:3,7:12,16:18];
Index0      = [4:6,13:15];
IndexSize   = 4:18;
for i = 1:length(files)
    if      ismember(i,Index9999)
        A = importdata(['/home/glaciology1/Documents/Data/GlacierTopos/', files(i).name],' ',6);
        A = double(A);  A(A==-9999) = NaN; 
    elseif  ismember(i,Index0)   
        A = importdata(['/home/glaciology1/Documents/Data/GlacierTopos/', files(i).name],' ',8);
        A = double(A);  A(A==0) = NaN; 
    end
    
    A(~any(~isnan(A), 2),:) =[];
    A(:,~any(~isnan(A), 1)) =[];
    if      ismember(i,IndexSize)
        A = A(2:end,:);
    end
    v(i,1) = cellstr(files(i).name(1:end-4)); eval([v{i,1} '= A;']);
end

%Create structure with full raster of data
for i = 1:length(v)/3
    param = char(v(i*3)); param = param(1:end-3);
    eval(['topo_full.G4.(param) = ',v{3*i-1,1},';']);
    eval(['topo_full.G2.(param) = ',v{3*i-2,1}],';');
    eval(['topo_full.G13.(param) = ',v{3*i,1}],';');
end
    
% %Sx is one row too big. Remove the first one -> doesn't exactly match the other params...
   topo_full.G4.Sx      = topo_full.G4.Sx(2:end,:);
   topo_full.G2.Sx      = topo_full.G2.Sx(2:end-1,:);
   topo_full.G13.Sx     = topo_full.G13.Sx(2:end,:);
   
        clear aspect* elev* north* profil* slope* Sx* tangent* files A i param v glacier ans Index*

%% Import Sampled Topo Params

%Import topos
[topo, ~, ~] = xlsread('SPOT_TopoParams_noZZ.xlsx','Sheet1','A1:F2353');

%Remove zigzag values (they skew the MLR results)
run OPTIONS.m; options.ZZ = 2;
run MAIN.m

%Get divisions for each glacier of the continous file imported
div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];

%Create structure of topo params for each glacier
for i = 1:3
name = char(options.glacier(i));
    slope            = topo(div(i,1):div(i,2),2);
    profileCurve     = topo(div(i,1):div(i,2),3);
    tangentCurve     = topo(div(i,1):div(i,2),4);
    aspect           = topo(div(i,1):div(i,2),5);
    elevation        = topo(div(i,1):div(i,2),6);
    
    topo_sampled.(name) = struct('aspect',aspect, 'elevation',elevation,...
        'profileCurve',profileCurve, 'slope',slope,'tangentCurve',tangentCurve);
end

% Sx import
[d300, d300text] = xlsread('d300_h0.xlsx','sampling_d300_h0','B1:BU2353');
[d200, d200text] = xlsread('d200_h0.xlsx','sampling_d200_h0','B1:BU2353');
[d100, d100text] = xlsread('d100_h0.xlsx','sampling_d100_h0','B1:BU2353');
    d300text = strcat('d300',d300text(1,:));
    d200text = strcat('d200',d200text(1,:));
    d100text = strcat('d100',d100text(1,:));
    
    %Sx correlation to SWE
    for i = 1:3
        name    = char(options.glacier(i));
        X       = [d100(div(i,1):div(i,2),:), d200(div(i,1):div(i,2),:), d300(div(i,1):div(i,2),:)]; %get Sx for all distances
        y       = SWE(i).swe; %get swe data
        text    = [d100text, d200text, d300text];
        
        CC      = corr([y,X]);    CC = CC(2:end,1);
        
        best    = find(abs(CC)==max(abs(CC)));  best = best(1,1);
      
        winddir(i,:)            = [text(best), num2cell(CC(best))]; %create structure with info on best wind direction and distance
        topo_sampled.(name).Sx  = X(:,best); %add Sx with best correlation to structure
    end
            clear best d300* d200* d100* i name text X* y CC

%Distance from centreline import
run CentrelineDistance.m

    clear aspect centreline corner distance dive elevation G glacier i profileCurve slope Sx tangentCurve topo X Y

%% Calculations

%Aspect correction -> 0degrees is N and goes clockwise
    for i = 1:3
        glacier = char(options.glacier(i));
        
        aspect = topo_full.(glacier).aspect;
            aspect = aspect -90;    aspect(aspect<0) = aspect(aspect<0)+360;   aspect = 360-aspect; 
        topo_full.(glacier).aspect = aspect;
        
        aspect = topo_sampled.(glacier).aspect;
            aspect = aspect -90;    aspect(aspect<0) = aspect(aspect<0)+360;   aspect = 360-aspect;             
        topo_sampled.(glacier).aspect = aspect;
    end

%Northness
    for i = 1:3
        glacier = char(options.glacier(i));
        topo_full.(glacier).northness    = cosd(topo_full.(glacier).aspect).*sind(topo_full.(glacier).slope);
        topo_sampled.(glacier).northness = cosd(topo_sampled.(glacier).aspect).*sind(topo_sampled.(glacier).slope);
    end

%Aspect -> get North South components
    for i = 1:3
        glacier = char(options.glacier(i));
        topo_full.(glacier).aspect    = cosd(topo_full.(glacier).aspect);
        topo_sampled.(glacier).aspect = cosd(topo_sampled.(glacier).aspect);
    end    
    
    clear aspect
%% Standardizing variables

%Keep a copy of non-standardized variables for range of params sampled
%plots
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

order = {'centreD','elevation','aspect','slope','northness','profileCurve',...
            'tangentCurve','Sx'};
for i = 1:3
   name     = char(options.glacier(i));
   topo_full.(name)         = orderfields(topo_full.(name),order);
   topo_full_ns.(name)      = orderfields(topo_full_ns.(name),order);   
   topo_sampled.(name)      = orderfields(topo_sampled.(name),order);
   topo_sampled_ns.(name)   = orderfields(topo_sampled_ns.(name),order);
end
    clear i name order
    
%% Get topo params for zigzag locations
%Put zigzags back in the final SWE structure (not related to topo)
% run OPTIONS.m; options.ZZ = 1;
% run MAIN.m
% 
% topo_sampled_wZZ = topo_sampled;
% 
% %Find the zigzag label in the sampled topo structure and match it with the 
% %zigzag values that are missing for the full matrix
% for i = 1:3;
%     name = char(options.glacier(i));
% 
%     zz = SWE(i).pattern =='ZZ';
%     zz_lab = char(SWE(i).label(:)); zz_lab = cellstr(zz_lab(:,1:8));
% 
%     zz_vals = find(~cellfun(@isempty,strfind(cellstr(SWE(i).label),'SWE')));
%     zz_valsName = char(SWE(i).label(zz_vals)); zz_valsName = zz_valsName(:,1:8);
% 
%     params = fieldnames(topo_sampled_wZZ.(name));
%     for k = 1:length(params)
%         topo = char(params(k));
%         
%         for j = 1:size(zz_valsName,1)
%             TT = ~cellfun(@isempty,strfind(zz_lab,zz_valsName(j,:)));
%             topo_sampled_wZZ.(name).(topo)(TT) = topo_sampled.(name).(topo)(zz_vals(j));
%         end
% 
%     end
% end
% 
%     clear j k fields params topo TT zz*


