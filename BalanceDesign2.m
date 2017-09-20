
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
 sizeR = 50;

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

%% Linear Regression

    namesP = fieldnames(pWB);
    nRuns = 100;
for p = 5%1:length(namesP)
    
    for mc = 19:nRuns;
    %Add some noise
    WBinput = WBnoise(pWB.(namesP{p}));
    
    %Linear regresion
    subsetWB.(namesP{p})(mc) = LinearRegression( WBinput, pTOPO.(namesP{p}), topo_full );
    end
    
    %Average WB distribution
        T = struct2table(subsetWB.(namesP{p}));
    for g = 1:3;    glacier = options.glacier{g};
        TT = T.(glacier);
        TT = reshape(TT,1,1,nRuns);
        TT = cell2mat(TT);
    subsetWBavg.(namesP{p}).(glacier) = mean(TT,3);
    end
end

for p = 1:length(namesP)
    %Plot resulting WB distribution from LR
    figure(p)
    PlotTopoParameter(subsetWBavg.(namesP{p}),'WB','WB (m w.e.)',pUTM.(namesP{p}), 'black', 'massB')
end

%% Difference Maps & RMSE

load Full.mat fullLR

    namesP = fieldnames(pWB); 
for p = 1:length(namesP)
    for g = 1:3;    glacier = options.glacier{g};
    %Difference map
    pDIFFmap.(namesP{p}).(glacier) = fullLR.S2.(glacier) - subsetWBavg.(namesP{p}).(glacier);
    
    %RMSE for each run
    for mc = 1:nRuns;
    diff = subsetWB.(namesP{p})(mc).(glacier)-fullLR.S2.(glacier);
    pRMSE.(glacier)(mc,p) = sqrt(nanmean((diff(:).^2)));
    end
        pRMSE.(glacier) = array2table(pRMSE.(glacier),'VariableNames',namesP);
    end
end

%Histogram of RMSE distributions
clf
for g = 1:3;    glacier = options.glacier{g};
    subplot(1,3,g)
for i = 5:-1:1
histogram(pRMSE.(glacier){:,i},0:0.004:0.22,'FaceAlpha',0.5, 'EdgeColor','none'); hold on
end
    legend(namesP(5:-1:1))
    xlabel('RMSE (m w.e.)'); ylabel('Frequency')
end

%% Get cell num raster

% [raster, info] = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/Donjek_7x7aspect.tif');
% 
% cellN = 1:size(raster,1)*size(raster,2);
% raster = reshape(cellN,[size(raster,1),size(raster,2)]);
% 
% geotiffwrite('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/CellNumber.tif',raster,info,'CoordRefSysCode','EPSG:32607')
%     clear
% % In QGIS: 1) convert sampling design to line, 2) convert to points, 3)
% % point sample the cell num raster, 4) Export to csv





