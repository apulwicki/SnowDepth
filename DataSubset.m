function [ SWEdata, TOPOdata ] = DataSubset( subset, clt, input )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global options
SWE             = input.SWE; 
topo_sampled    = input.topo_sampled; 
topo_sampled_ns = input.topo_sampled_ns;

for g = 1:3;
glacier = char(options.glacier(g));

% Pattern index %
if strcmp(subset, 'pattern')
    
     % Accumulation area
 if ischar(clt)     
     if strcmp(clt,'accum')
       AA.G4  = [150,364,579,644];
       AA.G2  = 753:756;
       AA.G13 = [922:929,931,932];
       I = false(size(SWE(g).swe));
       I(AA.(glacier)) = 1;
     end
    
 else   
    if clt == 1;     pattern = [{'UM'};{'LM'}];               %centreline
    elseif clt == 2; pattern = [{'UM'};{'LM'}];               %centreline + 4 transects
    elseif clt == 3; pattern = [{'UM'};{'LM'};{'UT'}];        %centreline + 3 transects
    elseif clt == 4; pattern = [{'UH'};{'LH'}];               %hourglass 
    elseif clt == 5; pattern = [{'UH'};{'LH'};{'UC'};{'LC'}]; %hourglass and circle
        %pattern = [{'LM'};{'LH'};{'LC'}];        %lower ablation
        %pattern = [{'UM'};{'UH'};{'UC'};{'UT'};{'BT'};]; %upper ablation
    end
    
    for p = 1:length(pattern)
        I1(:,p) = SWE(g).pattern == pattern(p,1); end
    I1 = any(I1,2);
    
     %Hourglass transects
    if clt == 2
            HGT = [83:92,111:119,21:36,49:56,...
                    240:248,267:275,371:378,450:458,407:413,483:488,...
                    589:602,629:644,745:760,785:792];
            HGT = categorical([HGT+0.1;HGT+0.2;HGT+0.3]);  HGT = HGT(:);
        for p = 1:length(HGT)
        I2(:,p) = SWE(g).label == HGT(p); end
        I2 = any(I2,2);
    I = any([I1, I2],2);
    elseif clt == 3
            HGT = [83:92,49:56,...
                    240:248,407:413,483:488,...
                    589:602,785:792];
            HGT = categorical([HGT+0.1;HGT+0.2;HGT+0.3]);  HGT = HGT(:);
        for p = 1:length(HGT)
        I2(:,p) = SWE(g).label == HGT(p); end
        I2 = any(I2,2);
    I = any([I1, I2],2);
    else I = I1;
    end
 end
 
% Sampling density %
elseif strcmp(subset, 'density')
    label = cellstr(SWE(g).label);
     %one/two/three people
    I1 = ~cellfun(@isempty,strfind(label,'.1'));
    I2 = ~cellfun(@isempty,strfind(label,'.2'));
    I3 = ~cellfun(@isempty,strfind(label,'.3'));
    if option.people == 1;      In = any(I1,2);
    elseif option.people == 2;  In = any([I1, I2],2);  
    elseif option.people == 3;  In = any([I1,I2,I3],2); 
    end
    
     %measurement location density
    T = strfind(label,'.');     T = cellfun(@isempty,T);
    label(T) = {'000.0'}; 
    label = char(label);
    if g == 1;
        T = label(:,3) == '.';      label(T,:) = strcat('0',label(T,:));
        T = label(:,2) == '.';      label(T,:) = strcat('00',label(T,:));
    end
    label = str2num(label);     label = floor(label);
    
    if option.density == 2;     wps = 4:2:max(label);
    elseif option.density == 3; wps = 4:3:max(label);
    end
    Id = ismember(label,wps);
    
     %combine
    I = In & Id;
    
elseif strcmp(subset, 'topoparam')
    
    if strcmp(option.lessgreat, 'less')
        I = topo_sampled_ns.(glacier).(option.topo) < option.value;
    elseif strcmp(option.lessgreat, 'greater')
        I = topo_sampled_ns.(glacier).(option.topo) > option.value;
    end
    
elseif strcmp(subset, 'random')
    
    n = option;
    I = randi(length(SWE(g).swe),n,1);
    
end

% Select SWE and topo data %
     %SWE data
    SWEdata.(glacier) = [SWE(g).swe(I), SWE(g).utm(I,1:2), SWE(g).cellN(I)];

     %Topo data
    param = fieldnames(topo_sampled.G4);
    for t = 1:length(param)
       P = char(param(t));
       TOPOdata.(glacier).(P) = topo_sampled.(glacier).(P)(I);
    end
    
clear label T I*
end 

end

