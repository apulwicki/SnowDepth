%% Universal Kriging with Monte Carlo Uncertainty Analysis
    %Alex Pulwicki, Jan 2018
%This script runs a Monte Carlo analysis of winter-balance data using the 
%DiceKriging package in R. 

%Source of uncertainty: 
%   sigma_GS = variability in gridcell-averaged b_w values, implemented by 
%               adding randomly chosen b_w from normal distribution to
%               input data
%   sigma_rho = get best estimate of distributed b_w for all density
%               options the take std of those
%   sigma_INT = obtain std of b_w estimate for each gridcell from the
%               DiceKriging package and find total std (sqrt of average
%               variance
%   sigma_ALL = mean and total std of data from MC


%% Generating input data for universal kriging

    %Loads input data: point-scale b_w, topographic parameters from sampled
    %cells, and options
    load TopoSWE.mat fullSWE topo_sampled
    run OPTIONS

    %Generates gridcell-averaged b_w values for input to the UK algorithm
    for d = 1:8;    den = options.DenOpt{d};
    [ inputSWE.(den), TOPOdata ] = ObsInCell( fullSWE.(den).input, topo_sampled);
    end

    %Generates normal distribution of b_w values based on zigzag std for each
    %glacier and a set of randomly chosen values from this distribution to add
    %to the input data (sigma_GS)
    varSIG = [0.027, 0.035, 0.04];      %zigzag std
    for g = 1:3;    glacier = options.glacier{g};
        varPD.(glacier)   = makedist('Normal','mu',0,'sigma',varSIG(g));    %normal distributions
        varSWE.(glacier)  = random(varPD.(glacier), length(inputSWE.S1.(glacier)),1000);   %set of b_w values to add 
    end

    %Initializing variables for uncertainty plotting (table format)   
    UKsigmaGS.mean = table(zeros(8,1),zeros(8,1),zeros(8,1),'VariableNames',options.glacier,'RowNames',options.DenOpt);
        UKsigmaGS.std  = UKsigmaGS.mean;      stdDen = UKsigmaGS.mean;
    UKsigmaINT.mean = UKsigmaGS.mean;     UKsigmaINT.std = UKsigmaGS.mean;
    UKsigmaALL     = table(zeros(2,1),zeros(2,1),zeros(2,1),'VariableNames',options.glacier,'RowNames',{'mean','std'});
    UKsigmaRHO     = UKsigmaALL;     

            clear g glacier varPD den d fullSWE topo_sampled

%% Basic kriging (no Monte Carlo) to get distributed b_w        
        
clc; format shortg; clock
    for d = 1:8;        den = options.DenOpt{d};
        display(den)    %Displays which density option the code is on at the moment 
    %Kriging
      fullUK.(den) =  KrigingR_G( inputSWE.(den) );
    end
clock

%Getting mean and std of Bw for plotting
%(***sigma_INT***)
    for d = 1:8;    den = options.DenOpt{d};
    for g = 1:3;    glacier = options.glacier{g};
        UKsigmaINT.mean{d,g}  = nanmean(fullUK.(den).(glacier).pred(:));
        %sigmaINT.std{d,g}   = sqrt(nanmean(fullUK.(den).(glacier).std(:).^2));  %glacier-wide std = sqrt of average variance  
    end
    end
    
%(***sigma_RHO***)
    for g = 1:3; glacier = options.glacier{g};
        UKsigmaRHO{1,g} = mean(UKsigmaINT.mean.(glacier));
        UKsigmaRHO{2,g} = std(UKsigmaINT.mean.(glacier));
    end
    clear d den g glacier

save('MonteCarloUKtemp.mat','-v7.3')
    
%% Kriging of the input data with Monte Carlo
%  # multistarts = 
%  universal kriging prediction with linear trend in easting and northing

numMC = 5;    %Number of Monte Carlo runs (paper says 1000)

clc; format shortg; clock
    for d = 3%1:8
        den = options.DenOpt{d};
        display(den)    %Displays which density option the code is on at the moment 

    for mc = 1:numMC
        if floor(mc/10)==mc/10; display(num2str(mc)); end %Displays which MC run the code is on

    %Adds the sigma_GS variability to input data
        for g = 1:3;        glacier = options.glacier{g};
            dataSWE.(den).(glacier)      = inputSWE.(den).(glacier);
            dataSWE.(den).(glacier)(:,1) = inputSWE.(den).(glacier)(:,1) + varSWE.(glacier)(:,mc);
            I = (dataSWE.(den).(glacier)(:,1)<0);
            dataSWE.(den).(glacier)(I,1) = 0;
        end

    %Kriging
      KRIGzz.(den)(mc) =  KrigingR_G( dataSWE.(den) );
    end
save('MonteCarloUKtemp.mat','KRIGzz','-v7.3') %save current data to temp file
clock
    end

%% Calculating glacier-wide winter balance for plotting purposes
    
for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};

%Calculate glacier-wide winter balance and std for all runs
%(***sigms_GS***)
    for mc = 1:numMC
        BwKRIG.(den).(glacier)(mc,1) = nanmean(KRIGzz.(den)(mc).(glacier).pred(:));
        BwKRIG.(den).(glacier)(mc,2) = sqrt(nanmean(KRIGzz.(den)(mc).(glacier).std(:).^2));
    end
    UKsigmaGS.mean{d,g} = mean(BwKRIG.(den).(glacier)(mc,1));
    UKsigmaGS.std{d,g}  = std(BwKRIG.(den).(glacier)(mc,1));
    
    stdDen{d,g} = mean(BwKRIG.(den).(glacier)(mc,2)); %std from interpolation and zigzags
end
%Bw and mean std for each density option
%(***sigms_ALL***)
    UKsigmaALL{1,g} = mean(UKsigmaGS.mean.(glacier)); 
    UKsigmaALL{2,g} = sqrt(mean(stdDen.(glacier).^2));  
end
    clear mc g d glacier den
    
%% Save final data set
save('MonteCarloUK.mat')