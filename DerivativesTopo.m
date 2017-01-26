   
    
%% Central Difference

h = 40; % step size
[f, R] = geotiffread('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Glacier Topo Maps/DonjekDEM_mini.tif');

%Smoothing grid
sizeG = 3; C = (sizeG+1)/2;

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


%% same cell

%  uiopen('/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/same_cell.csv',1) 

div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];

for g = 1:3
   glacier  = char(options.glacier(g));
   sameG = same_cell(div(g,1):div(g,2));

[A1, I]  = sort(sameG);

%sort everyone
swe     = SWE(g).swe(I);
topo    = struct2table(topo_sampled.(glacier));
topo    = topo(I,:);

T = diff(A1)==0;

   %not unique
   ind_no   = diff(A1)==0;
   swe_no   = swe(ind_no);
   sameG_no = unique(A1(ind_no));
   %unique
   ind_uni  = diff(A1)~=0;
   swe_uni  = swe(ind_uni);
   sameG_uni = A1(ind_uni);

for i = 1:length(sameG_no)
   ind      = sameG_no(i)==A1;
   swe_new(i,1) = mean(swe(ind));
   
   temp     = find(ind);   temp1    = temp(1,1);
   topo_new(i,:) = topo(temp1,:);
end
topo_new = topo(swe_new~=0,:);
swe_new  = swe_new(swe_new~=0,1);

swe = [swe_new;swe_uni];
topo = [topo_new;topo(ind_uni,:)];

data = [topo, table(swe,'VariableName',{'swe'})];
fitlm(data)
display(num2str(length(swe)))
end
   
   
   

J      = find(diff(A1)==0);

JJ=unique([J(:)',J(:)'+1])';

ss_r    = [false(1,1);ss(1:end-1)];
ss_mean = mean([swe(1:end-1),swe(2:end)],2);
swe(ss) = ss_mean(ss);
swe(ss_r) = [];

topo(ss_r,:) = [];
display(num2str(length(swe)))

data = [topo, table(swe,'VariableName',{'swe'})];

fitlm(data)

end



ss = same_cell(1:end-1)==same_cell(2:end);
ss_r = [false(1,1);ss(1:end-1)];

uniqueSWE = [SWE(1).swe; SWE(2).swe; SWE(3).swe];

    ss_mean = mean([uniqueSWE(1:end-1),uniqueSWE(2:end)],2);
    
    uniqueSWE(ss) = ss_mean(ss);
    uniqueSWE(ss_r) = [];


divN(1,:) = [1,             sum(~ss_r(div(1,1):div(1,2)))];
divN(2,:) = [divN(1,2)+1    divN(1,2)+sum(~ss_r(div(2,1):div(2,2)))];
divN(3,:) = [divN(2,2)+1    divN(2,2)+sum(~ss_r(div(3,1):div(3,2)-1))+1];


for g = 1:3
   glacier  = char(options.glacier(g));

   data = struct2table(topo_sampled.(glacier));
   data(ss_r(div(g,1):div(g,2)-1),:) = [];
   uniqueTOPO.(glacier) = [data,table(uniqueSWE(divN(g,1):divN(g,2)),'VariableName',{'swe'})];
        clear data
   fitlm(uniqueTOPO.(glacier))
   
end
clc