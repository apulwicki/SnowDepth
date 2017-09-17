
%% Get WB field that is the "true" field
    clear; close all
load Full.mat fullLR options
load TopoSWE.mat topo_full

 %Make matrix with utm of each grid cell
for g = 1:3;    glacier = options.glacier{g};
minE = min(options.rig.(glacier)(:,1));
minN = min(options.rig.(glacier)(:,2));

    nE = options.mapsize(g,2);            
    nN = options.mapsize(g,1);
    
    utmGridE.(glacier) = repmat([1:nE]*40+minE,nN,1);   
    utmGridN.(glacier) = repmat([nN:-1:1]'*40+minN,1,nE); 
end
 clear g* min* n*
%% Random choice

 % NUMBER OF RANDOM POINTS
 sizeR = 50;

for g = 1:3;    glacier = options.glacier{g};
data = fullLR.S2.(glacier);

 %Get random locations and the row,col for that cell
n = 1:length(data(:)); 
    n = reshape(n,size(data,1),size(data,2));
    n(isnan(data)) = NaN;
    n_NN = n(:);    n_NN(isnan(n_NN)) = [];
    
    I = randi(length(n_NN), sizeR, 1); 
Iloc = n_NN(I);

    row = nan(length(Iloc),1);  col=row;  WBsubval.(glacier)=row;
for i = 1:length(Iloc)
    [row(i), col(i)]=find(n==Iloc(i));
sampledWB.(glacier)(i,1) = data(row(i),col(i)); %Value of WB data at each random location
sampledUTM.(glacier)(i,2:3) = [utmGridE.(glacier)(row(i),col(i)), utmGridN.(glacier)(row(i),col(i))];
end
    clear n* I*

 %Get topoparam values at each random location
 param = fields(topo_full.G4);
for f = 1:length(param)
   field = param{f};
   for i = 1:length(row)
   sampledTOPO.(glacier).(field)(i,1) = topo_full.(glacier).(field)(row(i),col(i));
   end
end
end
 clear c* data f* g* i param row

%% Linear Regression

c = 1;  type = 'random';
subsetWB(c).(type) =  LinearRegression( sampledWB, sampledTOPO, topo_full );

PlotTopoParameter(subsetWB(c).(type),'WB','WB (m w.e.)',sampledUTM, 'black', 'massB')



%% Trying to get diff patterns
load Full.mat fullSWE fullLR
    clear sampled*
y = 5:1:70;
x = y+4;

for g = 1:3;    glacier = options.glacier{g};

    %Get UTM of grid cells sampled
    E = utmGridE.(glacier)(1,y);
    N = utmGridN.(glacier)(x,1);
sampledUTM.(glacier) = [ones(length(E),1),E',N];
    %Get WB values at each location
sampledWB.(glacier) = fullLR.S2.(glacier)(x,y);
    sampledWB.(glacier) = diag(sampledWB.(glacier));
    
     %Get topoparam values at each location
 param = fields(topo_full.G4);
for f = 1:length(param);   field = param{f};
   sampledTOPO.(glacier).(field) = topo_full.(glacier).(field)(x,y);
   sampledTOPO.(glacier).(field) = diag(sampledTOPO.(glacier).(field));
end

end

PlotTopoParameter(subsetWB(c).(type),'WB','WB (m w.e.)',sampledUTM, 'black', 'massB')

fitlm([struct2table(sampledTOPO.G13), table(sampledWB.G13)])

% c = 1;  type = 'centreline';
% subsetWB(c).(type) =  LinearRegression( sampledWB, sampledTOPO, topo_full );


%%

% - Use QGIS to get a set of utms for each pattern with a cell number 
% - Create structure with each pattern as a substructure with utm and cell
%number as entries
% - Create structure with WB at each pattern locations
pattern = 'LC';
clear e n wb
for g = 3;    glacier = options.glacier{g};
%Get UTM of measurements location within a pattern
    I = fullSWE.S2.(glacier).pattern == pattern;
    utm_pattern = fullSWE.S2.(glacier).utm(I,:);
        [~,ia,~] = unique(utm_pattern(:,3));
        ia = sort(ia);
    utm_pattern = utm_pattern(ia,:);
%Get UTM and WB for cells at the pattern locations
    e = zeros(length(utm_pattern),1); n = e; wb = e;
for i = 1:length(utm_pattern)
        [~, t1] =  min(abs(utmGridE.(glacier)(1,:) - utm_pattern(i,1)));
   e(i) = utmGridE.(glacier)(1,t1);
        [~, t2] =  min(abs(utmGridN.(glacier)(:,1) - utm_pattern(i,2))); 
   n(i) = utmGridN.(glacier)(t2,1); 
   wb(i) = fullLR.S2.(glacier)(t2,t1);
end

end

plot(n,e,'.')
%%



load TopoSWE.mat topo*

input.SWE = fullSWE.S2; input.topo_sampled = topo_sampled; input.topo_sampled_ns = topo_sampled_ns;
t = 1;
%     if     t == 1; type = 'centreline';          n = 10:5:55;    
%     elseif t == 2; type = 'CentreTransect4';     n = 10:10:100;  
%     elseif t == 3; type = 'CentreTransect3';     n = 10:10:100;
%     elseif t == 4; type = 'hourglass';           n = 10:10:100;  
%     elseif t == 5; type = 'hourglassCircle';     n = 10:10:100;  
%     elseif t == 6; type = 'circle';              n = 10:10:100;  
    

[ subsetSWE, TOPOdata ] = DataSubset( 'pattern', t, input );
[ subsetSWE, TOPOdata ] = ObsInCell(subsetSWE, TOPOdata); %=> utms of pattern points

%How do I find the grid cells in utmGrid that are closest to the subsetSWE
%utm coordinates?
T = find(utmGridE.G4 == subsetSWE.G4(1,2))
