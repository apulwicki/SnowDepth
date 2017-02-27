%% Selecting Data from Pattern

for g = 1:3;
glacier = char(options.glacier(g));
     %Pattern index**
%     pattern = [{'UM'};{'LM'}];               %centreline
%     pattern = [{'UH'};{'LH'};{'UC'};{'LC'}]; %hourglass and circle
%     pattern = [{'UT'};{'BT'}];               %transects
    pattern = [{'UM'};{'UH'};{'UC'};{'UT'};{'BT'};]; %upper ablation
%     pattern = [{'LM'};{'LH'};{'LC'}];        %lower ablation
    
    for p = 1:length(pattern)
        I.(glacier)(:,p) = SWE(g).pattern == pattern(p,1); end
    I.(glacier) = any(I.(glacier),2);

     %SWE data
    SWEdata.(glacier) = [SWE(g).swe(I.(glacier)), SWE(g).utm(I.(glacier),1:2)];

    %Topo data
    param = fieldnames(topo_sampled.G4);
    for t = 1:length(param)
       P = char(param(t));
       TOPOdata.(glacier).(P) = topo_sampled.(glacier).(P)(I.(glacier));
    end
end 

%% Plot - locations
%!!!*** fix map
    param = 'empty';
    topoParam.G4  = NaN(options.mapsize(1,:));
    topoParam.G2  = NaN(options.mapsize(2,:));
    topoParam.G13 = NaN(options.mapsize(3,:));
    topoParam.rig = options.rig;

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWEdata, 'colour')
    saveFIG('SamplingLocation_subset')

%% Regression Kriging
 
[ balanceRK, balance_residualsRK ] = RegressionKriging( SWEdata, TOPOdata, I, topo_full );

%% Plot estimate
 