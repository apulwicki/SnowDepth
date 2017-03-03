function [ sweRK, residualsRK, BMA ] = ...
    RegressionKriging( SWEdata, TOPOdata, topo_full, SWE )
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

%SWEdata = swe, easting, northing (all G)
%TOPOdata = topo_sampled but with selected data (all G)

global options
%% BMA
display('Working on BMA');

    cd BMS
    [BMAinit, BMAres] = BMS_R(SWEdata, TOPOdata);
    cd ..
    
for g = 1:3
    glacier = char(options.glacier(g));
BMA.(glacier) = BMAinit.(glacier);   
BMA.(glacier).Properties.VariableNames = {'BMACoefficient','BMAsemiR2','BMAunivarR2'};

residualsBMA.(glacier) = BMAres.(glacier);
end
        clear best i name X y t glacier stackMLR meanMLR g BMAinit BMAres t g glacier swe

% Predicting
   %check coeff order - BMS
    bmaCoeff = BMA.G4.Properties.RowNames(1:end-3);   topoCoeff = fieldnames(TOPOdata.G4);
    if ~isequal(bmaCoeff, topoCoeff)
        disp('Different order of coefficients between BMA and topo'); return; end
    
    for g = 1:3
    glacier = char(options.glacier(g));
         %Intercept
        sweBMA.(glacier) = repmat(BMA.(glacier){end-2,1}, size(topo_full.(glacier).centreD));
         %multiply coeffs and add them
        for n = 1:length(bmaCoeff)
            param               = char(bmaCoeff(n));
            sweT                = topo_full.(glacier).(param)*BMA.(glacier){n,1};
            sweBMA.(glacier) = sweBMA.(glacier) + sweT;
        end
         %Set min to 0
        sweBMA.(glacier)(sweBMA.(glacier)<0) = 0;
    end

%% Residuals Kriging - DiceKriging 
display('Working on kriging residuals');

for g = 1:3
        glacier = char(options.glacier(g));
    krig.(glacier) = KrigingR(residualsBMA.(glacier), SWEdata.(glacier)(:,2:3), glacier);
end    
    clear g glacier r
        
%% Combine
for g = 1:3
    glacier = char(options.glacier(g));
    sweRK.(glacier) = sweBMA.(glacier) + krig.(glacier).pred;
end   

 %Residuals of RK
for g = 1:3
    glacier = char(options.glacier(g));
    [~, I] = ismember(SWEdata.(glacier)(:,4), SWE(g).cellN);
T = sub2ind(size(sweRK.(glacier)),...
    floor(options.N.(glacier)(I)),floor(options.E.(glacier)(I)));

sampledRK.(glacier)      = sweRK.(glacier)(T); 
sampledBMA.(glacier)     = sweBMA.(glacier)(T);

residualsRK.(glacier)    = SWE(g).swe(I) - sampledRK.(glacier);
sweRK.LOO.(glacier)      = krig.(glacier).LOO + sampledBMA.(glacier);
end
    clear E g glacier N r T
    
display('Done')    
end

