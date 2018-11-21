%% Coefficients for all data 

% Selecting Data from Pattern

%     den = 'S2';
% 
% load TopoSWE.mat
% run OPTIONS
% clear DataObs_RMSE  
% 
% [ SWE, topo_sampled, STD ]    = ObsInCell(fullSWE.(den).input, topo_sampled);
% 
%     % Remove dc, aspect and Northness
% for g = 1:3;    glacier = options.glacier{g};
%     topo_sampled.(glacier) = rmfield(topo_sampled.(glacier),'centreD');
%     topo_sampled.(glacier) = rmfield(topo_sampled.(glacier),'aspect');
%     topo_sampled.(glacier) = rmfield(topo_sampled.(glacier),'northness');
% end
% 
% % BMA
% swe_input = SWE;
% topo_input = topo_sampled;
% 
% cd BMS
% [BMSinit, BMSres] = BMS_R(swe_input, topo_input);
% cd ..
% 
% for g = 1:3;        glacier = char(options.glacier(g));
% BMS.(glacier) = BMSinit.(glacier);   
% BMS.(glacier).Properties.VariableNames = {'BMSCoefficient','BMSsemiR2','BMSunivarR2'};
% residualsBMS.(glacier) = BMSres.(glacier);
% end
%% Coefficients for each pattern

% Selecting Data from Pattern

    den = 'S2';

load TopoSWE.mat
run OPTIONS

for t = [6,1,3,4,5,100]
if     t == 6; type = 'Circle';           subset = 'pattern';       
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

maxN = min([length(subsetSWE_temp.G4) length(subsetSWE_temp.G2) length(subsetSWE_temp.G13)]);

    % Remove dc, aspect and Northness
for g = 1:3;    glacier = options.glacier{g};
    TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'centreD');
    TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'aspect');
    TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'northness');
end



 for n = 5:maxN
     display([type, ' n=',num2str(n)])

for g = 1:3;    glacier = options.glacier{g};
    nI = randperm(maxN, n);
    %nI = floor(linspace(1,maxN,n));
    WBinput(n).(type).(glacier)   = subsetSWE_temp.(glacier)(nI,:);
        ff = fieldnames(TOPOdata_temp.(glacier));
    for i = 1:length(ff);    fname = ff{i};
    TOPOinput(n).(type).(glacier).(fname)  = TOPOdata_temp.(glacier).(fname)(nI,:); end

    % Correct the centreline values for invertable matrix when only centreline
    if g==3 && strcmp(type,'Centreline'); TOPOinput(n).(type).G13 =  rmfield(TOPOinput(n).(type).G13, 'centreD'); end
end

% Linear regression (BMA)
    swe_input = WBinput(n).(type);
    topo_input = TOPOinput(n).(type);
    
    cd BMS
    [BMSinit, BMSres] = BMS_R(swe_input, topo_input);
    cd ..
    
    for gg = 1:3;        glacier = char(options.glacier(gg));
    BMS(n).(type).(glacier) = BMSinit.(glacier);   
    BMS(n).(type).(glacier).Properties.VariableNames = {['BMSCoefficient_', num2str(t)],...
                                                 ['BMSsemiR2_', num2str(t)],...
                                                 ['BMSunivarR2_', num2str(t)]};
    residualsBMS(n).(type).(glacier) = BMSres.(glacier);    
    end

 end
end



%% Predict
% 
% for g = 1:3;    glacier = options.glacier{g};
%     sweMLR.(glacier) = repmat(MLR.(glacier)(1), options.mapsize(g,:));
%     mlrCoeff = MLR.(glacier)(2:end);    topoCoeff = fieldnames(topo_full.G4);
%         %multiply coeffs and add them
%     for m = 1:length(mlrCoeff)
%         param               = topoCoeff{m};
%         sweT                = topo_full.(glacier).(param)*mlrCoeff(m);
%         sweMLR.(glacier)    = sweMLR.(glacier) + sweT;
%     end
%         %Set min to 0
%     sweMLR.(glacier)(sweMLR.(glacier)<0) = 0;
% 
%     DataObs.(type).(glacier)(n,mc) = nanmean(sweMLR.(glacier)(:));
% 
%     RMSE
%     sampledtemp = sweMLR.(glacier)(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
%     estGrid     = diag(sampledtemp);
% 
%         %DataObs_RMSEeven.(type).(glacier)(n,mc) = sqrt(mean((estGrid-realGrid.(glacier)(:,1)).^2));
% end        