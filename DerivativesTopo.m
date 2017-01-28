   
    
%% Central Difference

h = 40; % step size
[f, R] = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/DonjekDEM_mini.tif');

for sizeG = 3:2:9
%Smoothing grid
%sizeG = 9; 
C = (sizeG+1)/2;

grid = ones(sizeG);
grid = grid/(sum(grid(:)));

f_filt = nan(size(f));
for i = C:size(f,1)-C
   for j =  C:size(f,2)-C
       T    = grid.*f(i-(C-1):i+(C-1),j-(C-1):j+(C-1));
       f_filt(i,j) = nansum(T(:));
   end
end
f_filt(f_filt<min(f(:))) = NaN;

%Derivative
    Fmn = f_filt(1:end-2,:);  Fme = f_filt(:,1:end-2);
    F0n = f_filt(2:end-1,:);  F0e = f_filt(:,2:end-1);  
    Fpn = f_filt(3:end,:);    Fpe = f_filt(:,3:end);
Y3n_1 = (Fpn-Fmn)/(2*h);            Y3e_1 = (Fpe-Fme)/(2*h);
Y3n_2 = (Fmn-2*F0n+Fpn)/h^2;        Y3e_2 = (Fme-2*F0e+Fpe)/h^2;
    Y3n_1 = [nan(1,size(Y3n_1,2));  Y3n_1; nan(1,size(Y3n_1,2))];
    Y3e_1 = [nan(size(Y3e_1,1),1),  Y3e_1, nan(size(Y3e_1,1),1)];
    Y3n_2 = [nan(1,size(Y3n_2,2));  Y3n_2; nan(1,size(Y3n_2,2))];
    Y3e_2 = [nan(size(Y3e_2,1),1),  Y3e_2, nan(size(Y3e_2,1),1)];
meanM  = (Y3n_1+Y3e_1)/2;
meanNE = (Y3n_2+Y3e_2)/2; 

meanM2(:,:,1) = abs(Y3n_1); meanM2(:,:,2) = abs(Y3e_1); 
meanM2(:,:,1) = nanmax(meanM2,[],3); meanM2(:,:,2) = [];

figure(2); clf
    subplot(2,2,1)
    imagesc(f); colorbar
        title('Original')
    
    subplot(2,2,2)
        imagesc(f_filt); colorbar; caxis([1400 3100])
        title(['Smoothed, window size = ',num2str(sizeG)]);
    
    subplot(2,2,3)
        imagesc(meanM2); colorbar; caxis([-1 1.5])
        title('1st Derivative Max')
    
    subplot(2,2,4)
        imagesc(meanNE); colorbar; caxis([-0.01 0.01])
        title('2nd Derivative Average')
 g = num2str(sizeG);
 saveFIG(['TopoDerivative',g,'x',g]);

%Save to geotiff 
location = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/Donjek_';   

filename_og     = [location,g,'x',g,'smooth.tif'];
filename_1n     = [location,g,'x',g,'slopeN.tif'];
filename_1e     = [location,g,'x',g,'slopeE.tif'];
filename_1mean  = [location,g,'x',g,'slopeMAX.tif'];
filename_2n     = [location,g,'x',g,'curveN.tif'];
filename_2e     = [location,g,'x',g,'curveE.tif'];
filename_2mean  = [location,g,'x',g,'curveMEAN.tif'];

info = geotiffinfo('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/DonjekDEM_mini.tif');

geotiffwrite(filename_og,f_filt,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_1n,Y3n_1,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_1e,Y3e_1,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_1mean,meanM2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_2n,Y3n_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_2e,Y3e_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_2mean,meanNE,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

end    
n = 1;
cell_num = zeros(size(f));
for i = 1:size(f,1)
    for j = 1:size(f,2)
    cell_num(i,j) = n;
    n = n+1;
    end
end

filename_cellnum = [location,'cell_num.tif'];
geotiffwrite(filename_cellnum,cell_num,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

 %% Finding highest correlation
 
%  uiopen('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/curve_sampled.csv',1) 
%  uiopen('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/slope_sampled.csv',1) 

 s = [1, length(SWE(1).swe)];  s(2,:) = s(1,2) + [1, length(SWE(2).swe)];
     s(3,:) = s(2,2) + [1, length(SWE(3).swe)];
     
swe =  [SWE(1).swe; SWE(2).swe;  SWE(3).swe];   

curvesampled(:,1:8) = [];   curvesampled(:,5:12) = [];
slopesampled(:,5:12) = [];   


allC    = corr(swe, curvesampled{:,:})';
allS    = corr(swe, slopesampled{:,:})';
G4C     = corr(SWE(1).swe, curvesampled{s(1,1):s(1,2),:})';
G2C     = corr(SWE(2).swe, curvesampled{s(2,1):s(2,2),:})';
G13C    = corr(SWE(3).swe, curvesampled{s(3,1):s(3,2),:})';
G4S     = corr(SWE(1).swe, slopesampled{s(1,1):s(1,2),:})';
G2S     = corr(SWE(2).swe, slopesampled{s(2,1):s(2,2),:})';
G13S    = corr(SWE(3).swe, slopesampled{s(3,1):s(3,2),:})';

curve_corr = table(allC, G4C, G2C, G13C, 'RowNames',curvesampled.Properties.VariableNames);
slope_corr = table(allS, G4S, G2S, G13S, 'RowNames',slopesampled.Properties.VariableNames);

    clear all* G*
 %%
clf    
figure(1)
% subplot(1,4,1)
%         D = [curvesampled{:,7},swe];
%     plot(D(:,1),D(:,2),'.'); hold on
%         T = fitlm(D(:,1),D(:,2));
% %     plot(D(:,1), T.Coefficients{1,1}+T.Coefficients{2,1}*D(:,2))
%     xlabel(char(curvesampled.Properties.VariableNames(7)))
%     ylabel('All SWE')
%     T.Rsquared.Ordinary
% subplot(1,4,2)
        D = [curvesampled{s(1,1):s(1,2),7},SWE(1).swe];
    plot(D(:,1),D(:,2),'.'); hold on
        T = fitlm(D(:,2),D(:,1));
%     plot(D(:,1), T.Coefficients{1,1}+T.Coefficients{2,1}*D(:,2))
    xlabel(char(curvesampled.Properties.VariableNames(7)))
    ylabel('G4 SWE')
    T.Rsquared.Ordinary
% subplot(1,4,3)
        D = [curvesampled{s(2,1):s(2,2),7},SWE(2).swe];
    plot(D(:,1),D(:,2),'.'); hold on
        T = fitlm(D(:,1),D(:,2));
%     plot(D(:,1), T.Coefficients{1,1}+T.Coefficients{2,1}*D(:,2))
    xlabel(char(curvesampled.Properties.VariableNames(7)))
    ylabel('G2 SWE')
    T.Rsquared.Ordinary
% subplot(1,4,4)
        D = [curvesampled{s(3,1):s(3,2),7},SWE(3).swe];
    plot(D(:,1),D(:,2),'.'); hold on
        T = fitlm(D(:,1),D(:,2));
%     plot(D(:,1), T.Coefficients{1,1}+T.Coefficients{2,1}*D(:,2))
    xlabel(char(curvesampled.Properties.VariableNames(7)))
    ylabel('G2 SWE')   
    T.Rsquared.Ordinary
    
    legend('G4','G2','G13')
    
    %%
    clf
    figure(2)
    subplot(2,2,1); i = 5;
    plot(curvesampled{s(1,1):s(1,2),i},SWE(1).swe,'.'); hold on
    plot(curvesampled{s(2,1):s(2,2),i},SWE(2).swe,'.'); hold on
    plot(curvesampled{s(3,1):s(3,2),i},SWE(3).swe,'.'); hold on
        title(char(curvesampled.Properties.VariableNames(i)))
        legend('G4','G2','G13')
        xlim([-300 300])
    
    subplot(2,2,2); i = 6;
    plot(curvesampled{s(1,1):s(1,2),i},SWE(1).swe,'.'); hold on
    plot(curvesampled{s(2,1):s(2,2),i},SWE(2).swe,'.'); hold on
    plot(curvesampled{s(3,1):s(3,2),i},SWE(3).swe,'.'); hold on
        title(char(curvesampled.Properties.VariableNames(i)))
        legend('G4','G2','G13')
        xlim([-300 300])
        
    subplot(2,2,3); i = 7;
    plot(curvesampled{s(1,1):s(1,2),i},SWE(1).swe,'.'); hold on
    plot(curvesampled{s(2,1):s(2,2),i},SWE(2).swe,'.'); hold on
    plot(curvesampled{s(3,1):s(3,2),i},SWE(3).swe,'.'); hold on
        title(char(curvesampled.Properties.VariableNames(i)))
        legend('G4','G2','G13')
        xlim([-300 300])
        
    subplot(2,2,4); i = 8;
    plot(curvesampled{s(1,1):s(1,2),i},SWE(1).swe,'.'); hold on
    plot(curvesampled{s(2,1):s(2,2),i},SWE(2).swe,'.'); hold on
    plot(curvesampled{s(3,1):s(3,2),i},SWE(3).swe,'.'); hold on
        title(char(curvesampled.Properties.VariableNames(i)))
        xlim([-300 300])
        legend('G4','G2','G13')        
%% same cell

%  uiopen('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/same_cell.csv',1) 

div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];
std_cell = [];

for g = 1:3
   glacier  = char(options.glacier(g));
   sameG = same_cell(div(g,1):div(g,2));

[A1, I]  = sort(sameG);

%sort everyone
swe     = SWE(g).swe(I);
topo    = struct2table(topo_sampled.(glacier));
topo    = topo(I,:);

T = diff(A1)==0;    T1 = [0;T(1:end-1)]; T = any([T,T1],2);
sameG_not = unique(A1(T));
for i = length(sameG_not):-1:1
   ind                  = find(sameG_not(i)==A1);
   swe(ind(1,1),1)      = mean(swe(ind));
   std_cell(i,g)        = std(swe(ind));
   swe(ind(2:end,1))    = [];
  
   topo(ind(2:end,1),:) = [];
   A1(ind(2:end,1))     = [];
end


data = [topo, table(swe,'VariableName',{'swe'})];
fitlm(data)
end
std_cell(std_cell==0) = NaN;

figure(1)
    boxplot(std_cell,'Labels',{'G4','G2','G13'})
    ylabel({'Standard Deviation with',' one DEM cell (m w.e.)'})


   
   