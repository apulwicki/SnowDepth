% Coefficients for each pattern

% Selecting Data from Pattern


type = 'Circle';            
% type = 'Centreline';           
% type = 'CentreTransect';   
% type = 'Hourglass';          
% type = 'HourCircle';       
% type = 'RandomSafe';     


% file_path = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/';
file_path = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/';

sample_size_list = 6:2:55;
%% Get synthetic snow distributions
load PaperII_Syntheic.mat

% Beta coeff ranges
elev_range  = [0.05, 0.810];
sx_range    = [-0.294, 0.260];
curv_range  = [-0.077, 0.02];
slope_range = [-0.09, 0.277];
intercept   = [0.6205, 0.2627, 0.2354];

% Normal random beta values
num_syn_models = 100;

elev_beta   = normrnd(mean(elev_range), std(elev_range), [1,num_syn_models]);
sx_beta     = normrnd(mean(sx_range),   std(sx_range),   [1,num_syn_models]);
curv_beta   = normrnd(mean(curv_range), std(curv_range), [1,num_syn_models]);
slope_beta  = normrnd(mean(slope_range),std(slope_range),[1,num_syn_models]);

    % Remove dc, aspect and Northness
    for g = 1:3;    glacier = options.glacier{g};
    topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
    topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
    end
    
for g=1:3; glacier = options.glacier{g};
AblationArea.(glacier)(AblationArea.(glacier)==-0.1)=NaN;
AblationArea.(glacier)(~isnan(AblationArea.(glacier)))=1;
end


for m = 1:num_syn_models

for g = 1:3;    glacier = options.glacier{g};
    syn_model   = repmat(intercept(g), options.mapsize(g,:));
    betaCoeff   = [elev_beta(m), sx_beta(m), curv_beta(m), slope_beta(m)];   
    topoCoeff   = fieldnames(topo_full.G4);
    
    for c = 1:length(betaCoeff)
        param     = topoCoeff{c};
        sweT      = topo_full.(glacier).(param)*betaCoeff(c);
        syn_model = syn_model + sweT;
    end
        
    syn_model(syn_model<0) = 0; %Set min to 0
    syn_model = syn_model.*AblationArea.(glacier);
    
    snowdist_model(m).(glacier) = syn_model;
    
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SAMPLING THEORETICAL FIELD FROM PATTERNS

% Get WB field that is the "true" field
   % clear; close all
% load TopoSWE.mat topo_full

 %Make matrix with utm of each grid cell
for g = 1:3;    glacier = options.glacier{g};
minE = min(options.rig.(glacier)(:,1));
minN = min(options.rig.(glacier)(:,2));

    nE = options.mapsize(g,2);            
    nN = options.mapsize(g,1);
    
    utmGridE.(glacier) = repmat([1:nE]*40+minE,nN,1);   
    utmGridN.(glacier) = repmat([nN:-1:1]'*40+minN,1,nE); 
end
%  clear g* min* n*

for mc = 1:num_syn_models

%Additing noise to distributed WB
fullWB = snowdist_model(mc);
 
% Obtaining pattern data for utm, topo, and wb

 %Get csv files
pattern.Circle = csvread([file_path,'CellNum_Circle.csv'],1,0);
pattern.Centreline = csvread([file_path,'CellNum_Centreline.csv'],1,0);
pattern.CentreTransect = csvread([file_path,'CellNum_Transverse.csv'],1,0);
pattern.Hourglass = csvread([file_path,'CellNum_Hourglass.csv'],1,0);
pattern.HourCircle = csvread([file_path,'CellNum_HourglassCircle.csv'],1,0);

namesP = {type}; %fieldnames(pattern);
for p = 1:length(namesP)
        [~,ia,~] = unique(pattern.(namesP{p})(:,3));
        ia = sort(ia);
    pattern.(namesP{p}) = pattern.(namesP{p})(ia,:);
    
 clear I*
for g = 1:3;    glacier = options.glacier{g};
        I_E(:,1) = pattern.(namesP{p})(:,1) > min(utmGridE.(glacier)(:));
        I_E(:,2) = pattern.(namesP{p})(:,1) < max(utmGridE.(glacier)(:));
    I_E = all(I_E,2);
        I_N(:,1) = pattern.(namesP{p})(:,2) > min(utmGridN.(glacier)(:));
        I_N(:,2) = pattern.(namesP{p})(:,2) < max(utmGridN.(glacier)(:));
    I_N = all(I_N,2);
    I = all([I_E I_N],2);
    
    utmPattern.(namesP{p}).(glacier)(:,1:3) = pattern.(namesP{p})(I,:); 
end

 %Get UTM and WB for cells at the pattern locations
 clear e n wb
 for g = 1:3;    glacier = options.glacier{g};
    e = zeros(length(utmPattern.(namesP{p}).(glacier)),1); num_models = e; wb = e;

    for i = 1:length(utmPattern.(namesP{p}).(glacier))
        [~, t1] =  min(abs(utmGridE.(glacier)(1,:) - utmPattern.(namesP{p}).(glacier)(i,1)));
       pUTM.(namesP{p}).(glacier)(i,1) = utmGridE.(glacier)(1,t1);
            [~, t2] =  min(abs(utmGridN.(glacier)(:,1) - utmPattern.(namesP{p}).(glacier)(i,2))); 
       pUTM.(namesP{p}).(glacier)(i,2) = utmGridN.(glacier)(t2,1); 
       
       pWB.(namesP{p}).(glacier)(i,1) = fullWB.(glacier)(t2,t1);
       
            topoparams = fieldnames(topo_full.G4);
       for f = 1:length(topoparams)
       pTOPO.(namesP{p}).(glacier).(topoparams{f})(i,1) = topo_full.(glacier).(topoparams{f})(t2,t1); 
       end
    end
end
end

% RANDOM - SAFE AREA   
sizeR = 200;

cells = dlmread([file_path,'SafeArea.csv']);

cells = cells(1:5:end-5,2:3);   
Gsplit = [1, 1950, 3235, length(cells)]; Gsplit = flip(Gsplit);
    for g = 1:3;    glacier = options.glacier{g};
        safeR.(glacier) = cells(Gsplit(g+1):Gsplit(g)-1,:);
            I = randi(length(safeR.(glacier)), sizeR, 1); 
        safeR.(glacier) = safeR.(glacier)(I,:);

    end

 %Get UTM and WB for cells at the pattern locations
 clear e n wb
 for g = 1:3;    glacier = options.glacier{g};
    e = zeros(length(safeR.(glacier)),1); num_models = e; wb = e;

    for i = 1:length(safeR.(glacier))
        [~, t1] =  min(abs(utmGridE.(glacier)(1,:) - safeR.(glacier)(i,1)));
       pUTM.RandomSafe.(glacier)(i,1) = utmGridE.(glacier)(1,t1);
            [~, t2] =  min(abs(utmGridN.(glacier)(:,1) - safeR.(glacier)(i,2))); 
       pUTM.RandomSafe.(glacier)(i,2) = utmGridN.(glacier)(t2,1); 
       
       pWB.RandomSafe.(glacier)(i,1) = fullWB.(glacier)(t2,t1);
       
            topoparams = fieldnames(topo_full.G4);
       for f = 1:length(topoparams)
       pTOPO.RandomSafe.(glacier).(topoparams{f})(i,1) = topo_full.(glacier).(topoparams{f})(t2,t1); 
       end
    end
 end

    
%%
namesP = {type}; %fieldnames(pWB);
%  namesP = {'Circle'};

real_measure = SampledCell(snowdist_model(mc));

for p = 1:length(namesP)

for ss = sample_size_list
   display([' Sample size: ',num2str(ss),' Pattern: ',namesP{p}, ' Run:',num2str(mc)])

   [WBinput, TOPOinput, UTMinput] = SubsetSampleSize( pWB, pTOPO, pUTM, ss );


    %Add some noise
    WBinputN = WBnoise(WBinput.(namesP{p}),'low');
    
        for g = 1:3;        glacier = char(options.glacier(g));
        
%         % BASIC LR
%         swe	    = WBinputN.(glacier)(:,1);
%         Xt      = struct2array(TOPOinput.(namesP{p}).(glacier));
%         X       = [ones(length(Xt),1), Xt];
        
%           % Get coefficients
%         coeffs = regress(swe, X);
% %         synMLR(ss).(namesP{p})(mc).(glacier)   = coeffs;
 
        % Linear regression (BMA)
        cd BMS
        [BMSinit, BMSres] = BMS_R_CI(WBinputN, TOPOinput.(namesP{p}));
        cd ..

        for gg = 1:3;        glacier = char(options.glacier(gg));
            coeffs = BMSinit.(glacier){:,1};
        BMS_CI.(glacier)(ss,gg,:) = coeffs;
               
        %Predict
        swe_pred = repmat(coeffs(5), options.mapsize(gg,:));
        betaCoeff = coeffs(1:4);    topoCoeff = fieldnames(topo_full.G4);
            %multiply coeffs and add them
        for num_models = 1:length(betaCoeff)
            param      = topoCoeff{num_models};
            sweT       = topo_full.(glacier).(param)*betaCoeff(num_models);
            swe_pred   = swe_pred + sweT;
        end
            %Set min to 0
        swe_pred(swe_pred<0) = 0;
        swe_pred = swe_pred.*AblationArea.(glacier); 
        
%         SynPred(ss).(namesP{p})(mc).(glacier) = swe_pred;

        % Calculate RMSE
        true_dist = snowdist_model(mc).(glacier);
        pred_dist = swe_pred;
        
        
        synRMSE_CI.(namesP{p}).(glacier)(ss,mc) = sqrt(nanmean(nanmean((pred_dist-true_dist).^2)));
        end
        
        end

end
end
end

save('PaperII_AA_BMS_CI.mat', 'synRMSE_CI', 'BMS_CI')