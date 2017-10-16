%% SAMPLING THEORETICAL FIELD FROM PATTERNS

% Get WB field that is the "true" field
    clear; close all
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
pattern.circle = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Circle.csv',1,0);
pattern.centreline = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Centreline.csv',1,0);
pattern.trans = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Transverse.csv',1,0);
pattern.hourglass = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_Hourglass.csv',1,0);
pattern.hourCircle = csvread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/CellNum_HourglassCircle.csv',1,0);

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
pWB.random.(glacier)(i,1) = data(row(i),col(i)); %Value of WB data at each random location
pUTM.random.(glacier)(i,2:3) = [utmGridE.(glacier)(row(i),col(i)), utmGridN.(glacier)(row(i),col(i))];
end
    clear n* I*

 %Get topoparam values at each random location
 param = fields(topo_full.G4);
for f = 1:length(param)
   field = param{f};
   for i = 1:length(row)
   pTOPO.random.(glacier).(field)(i,1) = topo_full.(glacier).(field)(row(i),col(i));
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
       pUTM.SafeRandom.(glacier)(i,1) = utmGridE.(glacier)(1,t1);
            [~, t2] =  min(abs(utmGridN.(glacier)(:,1) - safeR.(glacier)(i,2))); 
       pUTM.SafeRandom.(glacier)(i,2) = utmGridN.(glacier)(t2,1); 
       
       pWB.SafeRandom.(glacier)(i,1) = fullWB.(glacier)(t2,t1);
       
            topoparams = fieldnames(topo_full.G4);
       for f = 1:length(topoparams)
       pTOPO.SafeRandom.(glacier).(topoparams{f})(i,1) = topo_full.(glacier).(topoparams{f})(t2,t1); 
       end
    end
 end
 
    clear c* data f* g* i param row Gsplit safeR
    clear e f g* i* I* n* p pattern t1 t2 topoparams utm* wb WB*


    
%% THEORETICAL - All n for WB

load Full.mat fullLR
    namesP = fieldnames(pWB);
    %namesP = {'hourglass'};
    namesPfull = {'Circle','Centreline','Centre&Transverse','Hourglass',...
                  'Hourglass&Circle','Random', 'Safe Random'};
nRuns = 100;

for p = 1:length(namesP)
for g = 1:3; glacier = options.glacier{g};

for ss = 8:length(pWB.(namesP{p}).(glacier))

   [WBinput(ss), TOPOinput(ss), UTMinput(ss)] = SubsetSampleSize( pWB, pTOPO, pUTM, ss );

    
    display(['Sample size: ',num2str(ss),' Pattern: ',namesP{p}])
    
    for mc = 1:nRuns;
    %Add some noise
    WBinputN = WBnoise(WBinput(ss).(namesP{p}),'low');
    
    %Linear regresion
        swe	    = WBinputN.(glacier)(:,1);
        Xt      = struct2array(TOPOinput(ss).(namesP{p}).(glacier));
        X       = [ones(length(Xt),1), Xt];

        % Get coefficients
        MLR.(glacier)            = regress(swe, X);
        
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
        
        T.(namesP{p}).(glacier)(ss,mc) = nanmean(sweMLR.(glacier)(:));
    end
end
end
end
        C =[     0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330;...
        0.0588    0.3490    0.1216];

  % Figure  
load Patterns.mat T_*
T = T_high; 
     clf; n = 1;
for p = 1:length(namesP)
for g = 1:3; glacier = options.glacier{g};
   numPoints = 8:length(pWB.(namesP{p}).(glacier));

meanWB = mean(T.(namesP{p}).(glacier)(8:end,:),2);
stdWB  = std(T.(namesP{p}).(glacier)(8:end,:),[],2);
N10    = find((stdWB(12:end)./meanWB(12:end))<0.1,1); N10 = N10+11;
    subplot(length(namesP),3,n)
%plot(sampleSize,WBt.(glacier).(namesP{p}),'Color',P(p).Color); hold on
%errorbar(sampleSize,meanWB,stdWB); hold on

plot(numPoints,meanWB,'LineWidth',3,'Color',C(p,:)); hold on
    upper = meanWB + stdWB;
    lower = meanWB - stdWB;
fill([numPoints flip(numPoints)],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.3,'EdgeColor','none')

plot([min(numPoints) max(numPoints)],[nanmean(fullLR.S2.(glacier)(:)),nanmean(fullLR.S2.(glacier)(:))],'--k')

if ~isempty(N10)
plot([N10 N10],[0 1.2],':k','LineWidth',2')
end
      
    title([glacier,' ',namesPfull{p}])
    if g ==1; ylabel('WB (m w.e.)'); end
    if n>15; xlabel('Sample size'); end
ylim([0 1.2])
xlim([0 150])
n = n+1;
end
end 
    
    
    