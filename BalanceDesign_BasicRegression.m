%% SAMPLING THEORETICAL FIELD FROM PATTERNS

% Get WB field that is the "true" field
   % clear; close all
load Full.mat fullLR options
load TopoSWE.mat topo_full

 %Make matrix with utm of each grid cell
for g = 1:3;    glacier = options.glacier{g};
minE = min(options.rig.(glacier)(:,1));
minN = min(options.rig.(glacier)(:,2));

    nE = options.mapsize(g,2);            
    nN = options.mapsize(g,1);
    
    utmGridE.(glacier) = repmat([1:nE]*40+minE,nN,1);   
    utmGridN.(glacier) = repmat([nN:-1:1]'*40+minN,1,nE); 
end
 clear g* min* n*
 
%Additing noise to distributed WB
fullWB = fullLR.S2;
%     for g = 1:3;    glacier = options.glacier{g};
%     noise = normrnd( 0, options.zzstd(g), size(fullWB.(glacier),1),size(fullWB.(glacier),2));
%     fullWB.(glacier) = fullWB.(glacier) + noise;
%     end
 
% Obtaining pattern data for utm, topo, and wb

 %Get csv files
% pattern.Circle = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Circle.csv',1,0);
% pattern.Centreline = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Centreline.csv',1,0);
% pattern.CentreTransect = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Transverse.csv',1,0);
% pattern.Hourglass = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Hourglass.csv',1,0);
% pattern.HourCircle = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_HourglassCircle.csv',1,0);
pattern.Circle = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Circle.csv',1,0);
pattern.Centreline = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Centreline.csv',1,0);
pattern.CentreTransect = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Transverse.csv',1,0);
pattern.Hourglass = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Hourglass.csv',1,0);
pattern.HourCircle = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_HourglassCircle.csv',1,0);

namesP = fieldnames(pattern);
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
    e = zeros(length(utmPattern.(namesP{p}).(glacier)),1); n = e; wb = e;

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
% figure(p)
% subplot(1,3,g)
% plot(pUTM.(namesP{p}).(glacier)(:,1),pUTM.(namesP{p}).(glacier)(:,2),'.')
end
end

% RANDOM PATTERN

 %Number of random points
 sizeR = 200;

for g = 1:3;    glacier = options.glacier{g};
data = fullWB.(glacier);

 %Get random locations and the row,col for that cell
p = 1:length(data(:)); 
    p = reshape(p,size(data,1),size(data,2));
    p(isnan(data)) = NaN;
    n_NN = p(:);    n_NN(isnan(n_NN)) = [];
    
    I = randi(length(n_NN), sizeR, 1); 
Iloc = n_NN(I);

    row = nan(length(Iloc),1);  col=row;  WBsubval.(glacier)=row;
for i = 1:length(Iloc)
    [row(i), col(i)]=find(p==Iloc(i));
pWB.Random.(glacier)(i,1) = data(row(i),col(i)); %Value of WB data at each random location
pUTM.Random.(glacier)(i,2:3) = [utmGridE.(glacier)(row(i),col(i)), utmGridN.(glacier)(row(i),col(i))];
end
    clear n* I*

 %Get topoparam values at each random location
 param = fields(topo_full.G4);
for f = 1:length(param)
   field = param{f};
   for i = 1:length(row)
   pTOPO.Random.(glacier).(field)(i,1) = topo_full.(glacier).(field)(row(i),col(i));
   end
end
end
    
% RANDOM - SAFE AREA    
cells = dlmread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/SafeArea.csv');
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
    e = zeros(length(safeR.(glacier)),1); n = e; wb = e;

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
 
    clear c* data f* g* i param row Gsplit safeR
    clear e f g* i* I* n* p pattern t1 t2 topoparams utm* wb WB*


    
%% THEORETICAL - All n for WB
load Patterns.mat
load TopoSWE.mat topo_full
    namesP = fieldnames(pWB);
    %namesP = {'hourglass'};
nRuns = 100;
%for i = 1:100;



for p = 1:length(namesP)
for g = 1:3 
    glacier = options.glacier{g};

    %%%%%%%
%     pTOPO.(namesP{p}).(glacier) = rmfield(pTOPO.(namesP{p}).(glacier),'centreD');
%     pTOPO.(namesP{p}).(glacier) = rmfield(pTOPO.(namesP{p}).(glacier),'aspect');
%     pTOPO.(namesP{p}).(glacier) = rmfield(pTOPO.(namesP{p}).(glacier),'northness');
%     pTOPO.(namesP{p}).(glacier) = rmfield(pTOPO.(namesP{p}).(glacier),'slope');
%     pTOPO.(namesP{p}).(glacier) = rmfield(pTOPO.(namesP{p}).(glacier),'curvature');
%     pTOPO.(namesP{p}).(glacier) = rmfield(pTOPO.(namesP{p}).(glacier),'Sx');

for ss = 8:50%length(pWB.(namesP{p}).(glacier))

   [WBinput(ss), TOPOinput(ss), UTMinput(ss)] = SubsetSampleSize( pWB, pTOPO, pUTM, ss );

    
    display(['Glacier:',glacier,' Sample size: ',num2str(ss),' Pattern: ',namesP{p}])
    
    for mc = 1:nRuns
    %Add some noise
    WBinputN = WBnoise(WBinput(ss).(namesP{p}),'high');
    
    %Linear regresion
        swe	    = WBinputN.(glacier)(:,1);
        Xt      = struct2array(TOPOinput(ss).(namesP{p}).(glacier));
        X       = [ones(length(Xt),1), Xt];

        % Get coefficients
        MLR.(glacier)            = regress(swe, X);
%             ElevationSyn.(namesP{p}).(glacier)(ss,mc) = MLR.(glacier)(2);
%             CentreSyn.(namesP{p}).(glacier)(ss,mc) = MLR.(glacier)(3);
%             WindSyn.(namesP{p}).(glacier)(ss,mc) = MLR.(glacier)(8);
%             SlopeSyn.(namesP{p}).(glacier)(ss,mc) = MLR.(glacier)(5);
%             CurveSyn.(namesP{p}).(glacier)(ss,mc) = MLR.(glacier)(7);
        
        %Predict
        sweMLR.(glacier) = repmat(MLR.(glacier)(1), options.mapsize(g,:));
        mlrCoeff = MLR.(glacier)(2:end);    topoCoeff = fieldnames(topo_full.G4);
            %multiply coeffs and add them
        for n = 1:length(mlrCoeff)
            param               = topoCoeff{n};
            sweT                = topo_full.(glacier).(param)*mlrCoeff(n);
            sweMLR.(glacier)    = sweMLR.(glacier) + sweT;
        end
            %Set min to 0
        sweMLR.(glacier)(sweMLR.(glacier)<0) = 0;
        
        SynObs_High.(namesP{p}).(glacier)(ss,mc) = nanmean(sweMLR.(glacier)(:));
        
        if ss == 40;
        fullSynObs_High.(namesP{p})(ss).(glacier)(:,:,mc) = sweMLR.(glacier);
        end
    end
end
end
end

% basicM(i,6) = mean(SynObs_High.(namesP{p}).(glacier)(40,:));
% basicS(i,6) = std(SynObs_High.(namesP{p}).(glacier)(40,:));
% end


    
%% DATA - WB for all n

% Selecting Data from Pattern

    den = 'S2';
    nRuns = 30;


load TopoSWE.mat
run OPTIONS
clear DataObs_RMSE  

realGrid    = ObsInCell(fullSWE.(den).input, topo_sampled);

 
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

 for n = 6:3:maxN
     display([type, ' n=',num2str(n)])

for g = 1:3;    glacier = options.glacier{g};

    for mc = 1:nRuns

    nI = randperm(maxN, n);
    %nI = floor(linspace(1,maxN,n));
    WBinput(n).(type).(glacier)   = subsetSWE_temp.(glacier)(nI,:);
        ff = fieldnames(TOPOdata_temp.(glacier));
    for i = 1:length(ff);    fname = ff{i};
    TOPOinput(n).(type).(glacier).(fname)  = TOPOdata_temp.(glacier).(fname)(nI,:); end

    % Correct the centreline values for invertable matrix when only centreline
    if g==3 && strcmp(type,'Centreline'); TOPOinput(n).(type).G13 =  rmfield(TOPOinput(n).(type).G13, 'centreD'); end
      
% Linear regression    
        swe	    = WBinput(n).(type).(glacier)(:,1);
        Xt      = struct2array(TOPOinput(n).(type).(glacier));
        X       = [ones(length(Xt),1), Xt];

        % Get coefficients
        MLR(n).(type).(glacier)            = regress(swe, X);
%         ElevationReal.(type).(glacier)(n,mc) = MLR.(glacier)(2);

%         %Predict
%         sweMLR.(glacier) = repmat(MLR.(glacier)(1), options.mapsize(g,:));
%         mlrCoeff = MLR.(glacier)(2:end);    topoCoeff = fieldnames(topo_full.G4);
%             %multiply coeffs and add them
%         for m = 1:length(mlrCoeff)
%             param               = topoCoeff{m};
%             sweT                = topo_full.(glacier).(param)*mlrCoeff(m);
%             sweMLR.(glacier)    = sweMLR.(glacier) + sweT;
%         end
%             %Set min to 0
%         sweMLR.(glacier)(sweMLR.(glacier)<0) = 0;
%         
%         DataObs.(type).(glacier)(n,mc) = nanmean(sweMLR.(glacier)(:));
%         
        %RMSE
%         sampledtemp = sweMLR.(glacier)(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
%         estGrid     = diag(sampledtemp);

        %DataObs_RMSEeven.(type).(glacier)(n,mc) = sqrt(mean((estGrid-realGrid.(glacier)(:,1)).^2));
        
    end
        
end
% figure; plot(WBinput(n).(type).G4(:,2),WBinput(n).(type).G4(:,3),'.')
% title(type)
 end
end


%% PLOTTING - see Plot_PaperII

%% Distance travelled of random design

%%%%%%%%%%%%%%%%%%%%%%% All n %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(1); %clf; clc

for g = 1:3;    glacier = options.glacier{g};

E = pUTM.Random.(glacier)(2:end,2);      N = pUTM.Random.(glacier)(2:end,3);
Es = E(1);    Ns = N(1);

Dtravel = zeros(length(E),1);

for i = 1:length(E)
dist = sqrt((E-Es).^2 + (N-Ns).^2);

[D,I] = min(dist);      Dtravel(i+1) = Dtravel(i) + D;
Es = E(I);      Ns = N(I);
E(I) = [];      N(I) = [];
end

Dtravel = Dtravel/1000;
plot(Dtravel); hold on

display([glacier, ' ', num2str(Dtravel(40))])

end

legend(options.glacier, 'Location','northwest')
%%

%%%%%%%%%%%%%%%%%%%%%%% One n, diff random spots %%%%%%%%%%%%%%%%%%%%%%%%%
n = 40; clc
    Dtravel = zeros(n,nRuns);

for g = 1:3;    glacier = options.glacier{g};
for r = 1:nRuns
    Eruns = repmat(pUTM.Random.(glacier)(:,2),1,nRuns);
    Nruns = repmat(pUTM.Random.(glacier)(:,3),1,nRuns);
    ch = randperm(size(Eruns,1),40);
Eruns = Eruns(ch);  Nruns = Nruns(ch); 

E = Eruns(2:end);      N = Nruns(2:end);
Es = Eruns(1);         Ns = Nruns(1);

for i = 1:length(E)
dist = sqrt((E-Es).^2 + (N-Ns).^2);

[D,I] = min(dist);      Dtravel(i+1,r) = Dtravel(i,r) + D;
Es = E(I);      Ns = N(I);
E(I) = [];      N(I) = [];
end

%Dtravel = Dtravel/1000;

end
%display([glacier, ' ', num2str(min(Dtravel(n,:))/1000),' ',num2str(max(Dtravel(n,:))/1000)])
display([glacier, ' ', num2str(mean(Dtravel(n,:))/1000)])

end
plot(Dtravel);


%% Basic LR on all data

%Get basic LR coeffs
for g = 1:3
    glacier = char(options.glacier(g));
    X = [struct2table(topo_sampled.(glacier)), table(SWE(g).swe,'VariableNames',{'swe'})];
    basicLR.(glacier) = fitlm(X);
end

% Plot with full LR from paper
cmap = cbrewer('qual','Dark2',3);
figure(1);  clf
load Full.mat fullLR 

    %full LR coeffs
betas = zeros(7,3,8);
for d = 1:8; den = options.DenOpt{d};
    betas(:,:,d) = fullLR.(den).coeff{1:7,:}; 
end
fullLR_coeff = mean(betas,3);

title_list = {'G4','G2','G13'};
for g=1:3
    coeffs = [basicLR.(glacier).Coefficients{2:8,1}, fullLR_coeff(:,g)];
    subplot(1,3,g)
    B = bar(coeffs, 'EdgeColor','none');
    for i = 1:2; B(i).FaceColor = cmap(i,:); end
    legend('Basic','Full')
    title(title_list(g))
    set(gca,'xticklabel',options.topoVars)
end
  