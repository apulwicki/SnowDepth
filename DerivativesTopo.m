   
    
%% Central Difference

h = 40; % step size
[f, R] = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/DonjekDEM_mini.tif');

%Smoothing grid
sizeG = 9; C = (sizeG+1)/2;

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

figure(2); clf
    subplot(2,2,1)
    imagesc(f); colorbar
        title('Original')
    
    subplot(2,2,2)
        imagesc(f_filt); colorbar; caxis([1400 3100])
        title(['Smoothed, window size = ',num2str(sizeG)]);
    
    subplot(2,2,3)
        imagesc(meanM); colorbar; caxis([-1 1.5])
        title('1st Derivative Average')
    
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
filename_1mean  = [location,g,'x',g,'slopeMEAN.tif'];
filename_2n     = [location,g,'x',g,'curveN.tif'];
filename_2e     = [location,g,'x',g,'curveE.tif'];
filename_2mean  = [location,g,'x',g,'curveMEAN.tif'];

info = geotiffinfo('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/DonjekDEM_mini.tif');

geotiffwrite(filename_og,f_filt,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_1n,Y3n_1,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_1e,Y3e_1,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_1mean,meanM,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_2n,Y3n_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_2e,Y3e_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filename_2mean,meanNE,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

    
 %% Finding highest correlation
 
%  uiopen('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/curve_sampled.csv',1) 
%  uiopen('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/slope_sampled.csv',1) 

 s = [1, length(SWE(1).swe)];  s(2,:) = s(1,2) + [1, length(SWE(2).swe)];
     s(3,:) = s(2,2) + [1, length(SWE(3).swe)];
     
swe =  [SWE(1).swe; SWE(2).swe;  SWE(3).swe];   


allC    = corr(swe, curvesampled{:,:}).^2';
allS    = corr(swe, slopesampled{:,:}).^2';
G4C     = corr(SWE(1).swe, curvesampled{s(1,1):s(1,2),:}).^2';
G2C     = corr(SWE(2).swe, curvesampled{s(2,1):s(2,2),:}).^2';
G13C    = corr(SWE(3).swe, curvesampled{s(3,1):s(3,2),:}).^2';
G4S     = corr(SWE(1).swe, slopesampled{s(1,1):s(1,2),:}).^2';
G2S     = corr(SWE(2).swe, slopesampled{s(2,1):s(2,2),:}).^2';
G13S    = corr(SWE(3).swe, slopesampled{s(3,1):s(3,2),:}).^2';

curve_corr = table(allC, G4C, G2C, G13C, 'RowNames',curvesampled.Properties.VariableNames);
slope_corr = table(allS, G4S, G2S, G13S, 'RowNames',slopesampled.Properties.VariableNames);



 
 
 
 
 