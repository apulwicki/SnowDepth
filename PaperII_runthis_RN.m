% Coefficients for each pattern

% Selecting Data from Pattern

    den = 'S2';

load TopoSWE.mat
run OPTIONS

% t = 6; type = 'Circle';           subset = 'pattern';       
% t = 1; type = 'Centreline';       subset = 'pattern';     
% t = 3; type = 'CentreTransect';   subset = 'pattern';   
% t = 4; type = 'Hourglass';        subset = 'pattern';  
% t = 5; type = 'HourCircle';       subset = 'pattern';
t = 100; type = 'RandomSafe';     subset = 'random';



input.SWE = fullSWE.(den); input.topo_sampled = topo_sampled; 
input.topo_sampled_ns = topo_sampled_ns;

[ subsetSWE_temp, TOPOdata_temp ] = DataSubset( subset, t, input );

[ subsetSWE_temp, TOPOdata_temp ] = ObsInCell( subsetSWE_temp, TOPOdata_temp ); 

% maxN = min([length(subsetSWE_temp.G4) length(subsetSWE_temp.G2) length(subsetSWE_temp.G13)]);
maxN = 45;
numRand = 30;

    % Remove dc, aspect and Northness
for g = 1:3;    glacier = options.glacier{g};
    TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'centreD');
    TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'aspect');
    TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'northness');
end

BMS_RandomSafe.G4 = zeros(maxN, numRand,7);
BMS_RandomSafe.G2 = zeros(maxN, numRand,7);
BMS_RandomSafe.G13 = zeros(maxN, numRand,7);

for n = 6:3:maxN

 for x = 1:numRand  
     disp([type, ' n=',num2str(n), ' run=',num2str(x)])
     
for g = 1:3;    glacier = options.glacier{g};
    nI = randperm(maxN, n);
    %nI = floor(linspace(1,maxN,n));
    WBinput.(glacier)   = subsetSWE_temp.(glacier)(nI,:);
        ff = fieldnames(TOPOdata_temp.(glacier));
    for i = 1:length(ff);    fname = ff{i};
    TOPOinput.(glacier).(fname)  = TOPOdata_temp.(glacier).(fname)(nI,:); end

end

% Linear regression (BMA)
    cd BMS
    [BMSinit, BMSres] = BMS_R_RN(WBinput, TOPOinput);
    cd ..
    
    for gg = 1:3;        glacier = char(options.glacier(gg));
    BMS_RandomSafe.(glacier)(n,x,:) = BMSinit.(glacier){:,1};
    end

 end
end
 
save('BMS_RandomSafe.mat','BMS_RandomSafe')

