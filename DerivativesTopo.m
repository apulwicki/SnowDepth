
%2D derivative
h = 40; % step size
f = topo_full_ns.G2.elevation;

order = 1;
Yn = diff(f,order,1)/h;
Ye = diff(f,order,2)/h;

figure(1)
    subplot(1,3,1)
    imagesc(f); colorbar
    
    subplot(1,3,2)
    imagesc(Yn); colorbar
    
    subplot(1,3,3)
    imagesc(Ye); colorbar
    
    
%% Central Difference

h = 40; % step size
[f, R] = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/DonjekDEM_mini.tif');

%Smoothing grid
sizeG = 3; C = (sizeG+1)/2;
grid = ones(C);
for i = 1:C;
    for j = 1:C
        grid(i,j) = 1;%sqrt(((C-i)*h)^2+((C-j)*h)^2);
    end
end
grid = [grid, fliplr(grid(:,1:C-1))];     grid = [grid; flipud(grid(1:C-1,:))];
grid = 1./grid;     grid(C,C) = 1;
grid = grid/(sum(grid(:)));

f_filt = nan(size(f));
for i = C:size(f,1)-C
   for j =  C:size(f,2)-C
       T    = grid.*f(i-(C-1):i+(C-1),j-(C-1):j+(C-1));
       f_filt(i,j) = nansum(T(:));
   end
end
f_filt(f_filt<min(f(:))) = NaN;

%Stencil N = 3
    Fmn = f_filt(1:end-2,:);  Fme = f_filt(:,1:end-2);
    F0n = f_filt(2:end-1,:);  F0e = f_filt(:,2:end-1);  
    Fpn = f_filt(3:end,:);    Fpe = f_filt(:,3:end);
Y3n_1 = (Fpn-Fmn)/(2*h);            Y3e_1 = (Fpe-Fme)/(2*h);
Y3n_2 = (Fmn-2*F0n+Fpn)/h^2;        Y3e_2 = (Fme-2*F0e+Fpe)/h^2;
    Y3n_2 = [nan(1,size(Y3n_2,2));  Y3n_2; nan(1,size(Y3n_2,2))];
    Y3e_2 = [nan(size(Y3e_2,1),1),  Y3e_2, nan(size(Y3e_2,1),1)];

figure(2); clf
    subplot(3,2,1)
    imagesc(f); colorbar
        title('Original')
    
    subplot(3,2,2)
        imagesc(f_filt); colorbar; caxis([1900 3100])
        title(['Smoothed, window size = ',num2str(sizeG)]);
    
    subplot(3,2,3)
        imagesc(Y3n_1); colorbar; caxis([-1 1.5])
        title('1st Derivative - N')
    
    subplot(3,2,4)
        imagesc(Y3e_1); colorbar; caxis([-1 1.5])
        title('1st Derivative - E')

    subplot(3,2,5)
        imagesc(Y3n_2); colorbar; caxis([-0.01 0.01])
        title('2nd Derivative - N')
    
    subplot(3,2,6)
        imagesc(Y3e_2); colorbar; caxis([-0.01 0.01])
        title('2nd Derivative - E')
 
figure(1)
    imagesc((Y3n_2(:,3:end)+Y3e_2(3:end,:))/2); colorbar; caxis([-0.01 0.01])
    title({'2nd Derivative (average)','Inverse-distance weighted 3x3 grid'})
    filename = ['2ndDerivative_grid',num2str(sizeG)];
    print([options.path1, filename],'-dpng','-r0');

    g = num2str(sizeG);
filename = ['/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/Donjek_',g,'x',g,'smooth.tif'];
filenamen2 = ['/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/Donjek_',g,'x',g,'smooth_n2.tif'];
filenamee2 = ['/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/Donjek_',g,'x',g,'smooth_e2.tif'];

info = geotiffinfo('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/DonjekDEM_mini.tif');

geotiffwrite(filename,f_filt,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filenamen2,Y3n_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);
geotiffwrite(filenamee2,Y3e_2,R,'GeoKeyDirectoryTag',info.GeoTIFFTags.GeoKeyDirectoryTag);

%% Matlab raster analysis

[Z, R] = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/Smoothed5x5/Donjek5x5_latlong.tif');
[ASPECT, SLOPE, gradN, gradE] = gradientm(Z, R);
    
figure(3)
subplot(2,2,1)
    imagesc(cosd(ASPECT)); colorbar
    title('Aspect')
subplot(2,2,2)
    imagesc(SLOPE); colorbar
    title('Slope')
subplot(2,2,3)
    imagesc(gradN); colorbar; caxis([-10 0])
    title('gradN')
subplot(2,2,4)
    imagesc(gradE); colorbar; caxis([0 10])    
    title('gradE')
    
    
 %%
 
%  uiopen('/home/glaciology1/Documents/QGIS/Data/sampled_curvetest.csv',1) 
% clear topo*

clear curve
 s = [1, length(SWE(1).swe)];  s(2,:) = s(1,2) + [1, length(SWE(2).swe)];
     s(3,:) = s(2,2) + [1, length(SWE(3).swe)];

for g = 1:3 
    glacier = char(options.glacier(g));
        curve.(glacier).E  = e2(s(g,1):s(g,2));
%        curve.(glacier).N  = n2(s(g,1):s(g,2));
       curve.(glacier).meanEN = mean([e2(s(g,1):s(g,2)),n2(s(g,1):s(g,2))],2);

      curve.(glacier).P  = Profilecu(s(g,1):s(g,2));
      curve.(glacier).T  = Tangential(s(g,1):s(g,2));
     curve.(glacier).meanPT = mean([Profilecu(s(g,1):s(g,2)),Tangential(s(g,1):s(g,2))],2);
 
end
%     clear e2 n2 Profilecu Tangential s g glacier

%Standardizing variables
params = fieldnames(curve.G4);
for i = 1:3
name = char(options.glacier(i));
    for t = 1:length(params)
    field = char(params(t));
    
    curve.(name).(field) = (curve.(name).(field)-...
        mean(curve.(name).(field)))/std(curve.(name).(field));
    end
end    
    





t=1;
for i = 1:3
    y       = SWE(i).swe;
    glacier = char(options.glacier(i)); 
        display(['glacier = ',glacier]);
    X       = curve.(glacier);

    MLRcurve.(glacier) = MLRcalval(y, X);

    cd BMS
    [BMSinit, BMSres] = BMS_R(SWE, topo_sampled);
    cd ..
    BMS(t).(glacier)     = BMSinit.(glacier);   

end

 subplot(1,2,1)
 plot(Profilecu,swe_opt8,'.')
  subplot(1,2,2)
 plot(Tangential,swe_opt8,'.')
    
 
 P_LR = fitlm(Profilecu,swe_opt8);
 figure; plot(P_LR)
 
 
 X = [ones(length(e2),1), mean([e2,n2],2)];
 X = [ones(length(e2),1), mean([Profilecu,Tangential],2)];
 X = [ones(length(e2),1), e2];
 X = [ones(length(e2),1), n2];
 X = [ones(length(e2),1), Profilecu];
 X = [ones(length(e2),1), Tangential];
 

y = [SWE(1).swe;SWE(2).swe;SWE(3).swe];
[b,~,~,~,stats] = regress(y,X);
 stats
 
 %%
 j = 3;
         glacier = char(options.glacier(j)); 
        X       = topo_sampled.(glacier);

 M = struct2table(X);

 %Everyone
 xE = [ones(size(M,1),1), M{:,:}];
 y = SWE(j).swe;

 [b,~, ~, ~, stats] = regress(y,xE);
 y_hatE = sum(repmat(b',size(xE,1),1).*xE,2);
 SSE_E = sumsqr(y_hatE-y);
 
 for i = 1:size(M,2)
     M1=M{:,:}; M1(:,i)=[]; 
 %No one
 xN = [ones(size(M1,1),1), M1];

 [bN,~, ~, ~, stats] = regress(y,xN);
 y_hatN = sum(repmat(bN',size(xE,1),1).*xN,2);
 SSE_N = sumsqr(y-y_hatN);
 
 partR(i) = (SSE_N-SSE_E)/SSE_N;
 end
 
 partR'
 
 r02 = corr(y,
 
 
 
 
 