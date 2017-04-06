function [ SWEdata, TOPOdata ] = DataSubset( subset, option, input )
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
    if option.clt == 1;     pattern = [{'UM'};{'LM'}];               %centreline
    elseif option.clt == 2; pattern = [{'UH'};{'LH'};{'UC'};{'LC'}]; %hourglass and circle
    elseif option.clt == 3; pattern = [{'UT'};{'BT'}];               %transects
    elseif option.clt == 4; pattern = [{'UM'};{'UH'};{'UC'};{'UT'};{'BT'};]; %upper ablation
    elseif option.clt == 5; pattern = [{'LM'};{'LH'};{'LC'}];        %lower ablation
    end
    
    for p = 1:length(pattern)
        I(:,p) = SWE(g).pattern == pattern(p,1); end
    I = any(I,2);

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

