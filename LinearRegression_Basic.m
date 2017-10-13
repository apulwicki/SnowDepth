function [ sweMLR ] = LinearRegression_Basic( SWEdata, TOPOdata, topo_full )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global options 
for i = 1:3
        glacier = options.glacier{i}; 
        swe	    = SWEdata.(glacier)(:,1);
        Xt      = struct2array(TOPOdata.(glacier));
        X       = [ones(length(Xt),1), Xt];

        % Get coefficients
        MLR.(glacier)            = regress(swe, X);
        
        %Predict
        sweMLR.(glacier) = repmat(MLR.(glacier)(1), options.mapsize(i,:));
        mlrCoeff = MLR.(glacier)(2:end);    topoCoeff = fieldnames(topo_full.G4);
            %multiply coeffs and add them
        for n = 1:length(mlrCoeff)
            param               = topoCoeff{n};
            sweT                = topo_full.(glacier).(param)*mlrCoeff(n);
            sweMLR.(glacier)    = sweMLR.(glacier) + sweT;
        end
            %Set min to 0
        sweMLR.(glacier)(sweMLR.(glacier)<0) = 0;
end

sweMLR.coeff = [MLR.G4(:,1), MLR.G2(:,1), MLR.G13(:,1)];

end


