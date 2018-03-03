%% Ordinary Kriging with Monte Carlo Uncertainty Analysis
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


%% Generating input data for ordinary kriging

    %Loads input data: point-scale b_w, topographic parameters from sampled
    %cells, and options
    load TopoSWE.mat fullSWE topo_sampled
    run OPTIONS

    %Generates gridcell-averaged b_w values for input to the OK algorithm
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
    OKsigmaGS.mean = table(zeros(8,1),zeros(8,1),zeros(8,1),'VariableNames',options.glacier,'RowNames',options.DenOpt);
        OKsigmaGS.std  = OKsigmaGS.mean;      stdDen = OKsigmaGS.mean;
    OKsigmaINT.mean = OKsigmaGS.mean;     OKsigmaINT.std = OKsigmaGS.mean;
    OKsigmaALL     = table(zeros(2,1),zeros(2,1),zeros(2,1),'VariableNames',options.glacier,'RowNames',{'mean','std'});
    OKsigmaRHO     = OKsigmaALL;     

            clear g glacier varPD den d fullSWE topo_sampled

%% Basic kriging (no Monte Carlo) to get distributed b_w        
        
clc; format shortg; clock
    for d = 1:8        
        den = options.DenOpt{d};
        display(den)    %Displays which density option the code is on at the moment 
    %Kriging
      fullOK.(den) =  KrigingR_G( inputSWE.(den) );
clock
    end


%Getting mean and std of Bw for plotting
%(***sigma_INT***)
    for d = 1:8;    den = options.DenOpt{d};
    for g = 1:3;    glacier = options.glacier{g};
        OKsigmaINT.mean{d,g}  = nanmean(fullOK.(den).(glacier).pred(:));
        OKsigmaINT.std{d,g}   = sqrt(nanmean(fullOK.(den).(glacier).std(:).^2));  %glacier-wide std = sqrt of average variance  
        
            X = 0:0.01:1.5;
        BwKRIGinterp.(den).(glacier) = normpdf(X, OKsigmaINT.mean{d,g}, OKsigmaINT.std{d,g});
    end
    end
    
%(***sigma_RHO***)
    for g = 1:3; glacier = options.glacier{g};
        OKsigmaRHO{1,g} = mean(OKsigmaINT.mean.(glacier));
        OKsigmaRHO{2,g} = std(OKsigmaINT.mean.(glacier));
    end
    clear d den g glacier

save('MonteCarloOKtemp.mat','-v7.3')
    
%% Kriging of the input data with Monte Carlo
%  # multistarts = 50
%  ordinary kriging prediction

numMC = 1000;    %Number of Monte Carlo runs (paper says 1000)

clc; format shortg; clock
    for d = 1
        den = options.DenOpt{d};
        display(den)    %Displays which density option the code is on at the moment 

    for mc = 501:numMC
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
save('MonteCarloOKtemp.mat','-v7.3') %save current data to temp file
clock
    end

%% Calculating glacier-wide winter balance for plotting purposes

    averagestdGS = OKsigmaGS.mean;
for g = 1:3;    glacier = options.glacier{g};
for d = 1:8;    den = options.DenOpt{d};

%Calculate glacier-wide winter balance and std for all runs
%(***sigms_GS***)
    for mc = 1:numMC
        BwKRIGzz.(den).(glacier)(mc) = nanmean(KRIGzz.(den)(mc).(glacier).pred(:));
        BwKRIGzz_std.(den).(glacier)(mc) = sqrt(nanmean(KRIGzz.(den)(mc).(glacier).std(:).^2));
               
    end
    OKsigmaGS.mean{d,g} = mean(BwKRIGzz.(den).(glacier));
    OKsigmaGS.std{d,g}  = std(BwKRIGzz.(den).(glacier));
    
    averagestdGS{d,g} = sqrt(mean(BwKRIGzz_std.(den).(glacier).^2));
end
%Bw and mean std for each density option
%(***sigms_ALL***)
    OKsigmaALL{1,g} = mean(OKsigmaGS.mean{:,g}); 
    OKsigmaALL{2,g} = sqrt(mean(averagestdGS{:,g}.^2));  
    
            X = 0:0.01:1.5;
    BwKRIGall.(glacier) = normpdf(X, OKsigmaALL{1,g}, OKsigmaALL{2,g});

end
    clear mc g d glacier den

    
%% Relative uncertainty calculation
% Calculates the differences in estaimted WB for each gridcell for the
% first 100 runs of the OK and then normalizes it. This corresponds to the 
% relative uncertainity (i.e. how much the value varies at each gridcell)

for d = 1:8 
                den = options.DenOpt{d};
for g = 1:3 
                glacier = options.glacier{g};
clear F G A H U W X %clearing needed temporary variables

%Initializing variables
    s = size(KRIGzz.(den)(1).(glacier).pred); %size of glacier grid
    n = 1; p = 1; % counters
    runs = 20;   % num of OK runs over which to calculate relative uncertainity (~100 to keep matrices reasonable size)
%Create stack of OK runs
    for i = 1:runs
    F(n:n+s(1)-1,:) = KRIGzz.(den)(i).(glacier).pred; %stack along first matrix dimension
    G(:,p:p+s(2)-1) = KRIGzz.(den)(i).(glacier).pred; %stack along second matrix dimension
    n = n+s(1); %increase counters by size of glacier grid
    p = p+s(2);
    end
%Repeat stack along other matrix dimension (gets all possible combinations
%of gridcells)
    F = repmat(F,1,runs);
    G = repmat(G,runs,1);

%Subtract two matrices to get all differences
    H = F-G;
    H = tril(H);   %Remove repeat data
    H(H<=0) = NaN; %Select only positives (negatives are repeats)
    U = isnan(H);  %Will be used to remove data in subsequent matrices

%Reshaping H and U to be a 3D matrix (glacier size by runs)
        n=1;
    for i = 1:s(1):size(H,1)
    for j = 1:s(2):size(H,2)
        A(:,:,n) = H(i:i+s(1)-1,j:j+s(2)-1);
        W(:,:,n) = U(i:i+s(1)-1,j:j+s(2)-1);
        n=n+1;
    end
    end

%Find where W (reshaped U, which is all nan values) is all nan and remove
%data
    for i = size(A,3):-1:1
    X(i) = all(all(W(:,:,i))); %Indices for where W is all nan
    end
    A(:,:,X) = []; %Remove all nan data

%Sum all differences and make it relative    
    D = nansum(A,3); %Sum all differences

    DOK.(den).(glacier) = D/max(D(:)); %Calculate difference in relation to max diff for each glacier
    DOK.(den).(glacier)(options.mapNaN.(glacier)) = NaN; %Set glacier outlines
end
end
        clear F G A H U W X d D den g glacier i j n p runs s ans %clearing needed temporary variables
        
%% Save final data set
save('MonteCarloOK.mat')