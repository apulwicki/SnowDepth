function [ sweMLR, CI ] = LinearRegression( SWEdata, TOPOdata, topo_full )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global options 
for i = 1:3
        glacier = options.glacier{i}; 
        swe	    = SWEdata.(glacier)(:,1);
        X       = TOPOdata.(glacier);

        % Get coefficients
        [MLR.(glacier), ~, CI.(glacier)]            = MLRcalval(swe, X);
        MLR.(glacier).Properties.VariableNames(1)   = options.glacier(i);
        
        %Predict
        mlrCoeff = MLR.G4.Properties.RowNames(1:end-3);   topoCoeff = fieldnames(topo_full.G4);
        sweMLR.(glacier) = repmat(MLR.(glacier){end-2,1}, size(topo_full.(glacier).centreD));
            %multiply coeffs and add them
        for n = 1:length(mlrCoeff)
            param               = mlrCoeff{n};
            sweT                = topo_full.(glacier).(param)*MLR.(glacier){n,1};
            sweMLR.(glacier)    = sweMLR.(glacier) + sweT;
        end
            %Set min to 0
        sweMLR.(glacier)(sweMLR.(glacier)<0) = 0;
end

sweMLR.coeff = [MLR.G4(:,1), MLR.G2(:,1), MLR.G13(:,1)];

end

