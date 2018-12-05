% Coefficients for each pattern

% Selecting Data from Pattern

    den = 'S2';

load TopoSWE.mat
run OPTIONS

for t = [6,1,3,4,5]
if t == 6; type = 'Circle';           subset = 'pattern';       
elseif t == 1; type = 'Centreline';       subset = 'pattern';     
elseif t == 3; type = 'CentreTransect';   subset = 'pattern';   
elseif t == 4; type = 'Hourglass';        subset = 'pattern';  
elseif t == 5; type = 'HourCircle';       subset = 'pattern';
elseif t == 100; type = 'RandomSafe';     subset = 'random';
end


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


for n = 5:3:maxN
     display([type, ' n=',num2str(n)])

 for x = 1:numRand    
     
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
    [BMSinit, BMSres] = BMS_R(WBinput, TOPOinput);
    cd ..
    
    for gg = 1:3;        glacier = char(options.glacier(gg));
    BMS.(type).(glacier)(n,x,:) = BMSinit.(glacier){:,1};
    end

 end
end
 

end

save('PaperII_AblationArea_Dec4.mat')
