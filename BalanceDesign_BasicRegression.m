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
 sizeR = 100;

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
    clear c* data f* g* i param row
    clear e f g* i* I* n* p pattern t1 t2 topoparams utm* wb WB*
    
    
%% Linear Regression & Points selection

%load Full.mat fullLR
    namesP = fieldnames(pWB); 
nRuns = 100;

numPoints = zeros(1,1); numPoints(1) = 8;
for n = 2:5;
    numPoints(n,1) = numPoints(n-1)*2-1;
end


for n = 1:5;
ss = numPoints(n);

   [WBinput(ss), TOPOinput(ss), UTMinput(ss)] = SubsetSampleSize( pWB, pTOPO, pUTM, ss );

for p = 1:length(namesP)
    
    display(['Sample size: ',num2str(ss),' Pattern: ',namesP{p}])
    
    for mc = 1:nRuns;
    %Add some noise
    WBinputN = WBnoise(WBinput(ss).(namesP{p}),'high');
    
    %Linear regresion
    subsetWB(ss).(namesP{p})(mc) = LinearRegression_Basic(WBinputN,...
                                            TOPOinput(ss).(namesP{p}),...
                                            topo_full);
    end
    
save('PatternsTempBasicLR_HighNoise.mat','*input','subset*','-v7.3')    
end
end

%% WB vs n
load Full.mat fullLR
run OPTIONS.m

clear WBt
    namesP = fieldnames(subsetWB); 
    namesPfull = {'Circle','Centreline','Centre&Transverse','Hourglass','Hourglass&Circle','Random'};
sampleSize = 40:5:100;
nRuns = 100;

for g = 1:3; glacier = options.glacier{g};

for s = 1:length(numPoints)
ss = numPoints(s);
for p = 1:length(namesP)
    for mc = 1:nRuns;
WBt.(glacier).(namesP{p})(s,mc) = nanmean(subsetWB(ss).(namesP{p})(mc).(glacier)(:));
    end
    WBt.(glacier).(namesP{p})(WBt.(glacier).(namesP{p})==0) = NaN;
end
end
    
end

    C =[     0    0.4470    0.7410;...
        0.8500    0.3250    0.0980;...
        0.9290    0.6940    0.1250;...
        0.4940    0.1840    0.5560;...
        0.4660    0.6740    0.1880;...
        0.3010    0.7450    0.9330];

% Plot - WB vs n
clf; n = 1;
for p = 1:length(namesP)
for g = 1:3; glacier = options.glacier{g};

meanWB = mean(WBt.(glacier).(namesP{p}),2);
stdWB  = std(WBt.(glacier).(namesP{p}),[],2);
    subplot(length(namesP),3,n)
%plot(sampleSize,WBt.(glacier).(namesP{p}),'Color',P(p).Color); hold on
%errorbar(sampleSize,meanWB,stdWB); hold on

plot(numPoints,meanWB,'-d','LineWidth',3,'Color',C(p,:)); hold on
    upper = meanWB + stdWB;
    lower = meanWB - stdWB;
fill([numPoints' flip(numPoints')],[upper',flip(lower')],...
     C(p,:),'FaceAlpha',0.3,'EdgeColor','none')

plot([min(numPoints) max(numPoints)],[nanmean(fullLR.S2.(glacier)(:)),nanmean(fullLR.S2.(glacier)(:))],'--k')

    title([glacier,' ',namesPfull{p}])
    if g ==1; ylabel('WB (m w.e.)'); end
    if n>15; xlabel('Sample size'); end
ylim([0 1.5])
n = n+1;
end
end
%legend(P, namesP)
    
    
    
    
    
    
    
    
    
    
    