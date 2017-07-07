%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Tutorial 4 (1 Oct 2015)
% PCA on real data
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Example 2:
% Gridded monthly sea surface temperature data for Tropical Pacific
% from ERA Interim reanalysis
% Period: Jan 1979 to Jun 2015
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

months={'Oct' 'Nov' 'Dec' 'Jan' 'Feb' 'Mar' 'Apr' 'May'};

% this is nc formated file (binary, structured file)
filename='VIWT_Oct-May_LowRes.nc';

% the following command will display the file content
ncdisp(filename);

% these are the factors read from the file
scale_factor  = 0.0012045;
add_offset    = 39.5782;
missing_value = -32767;

% opening the file
ncid = netcdf.open(filename,'NC_NOWRITE');   % this is the matlab function for opening nc files
[numdims, numvars, numglobalatts, unlimdimID] = netcdf.inq(ncid);  % this is all the structure in the nc file

% this loop looks for a particular variable from the set of variables 
for jj=1:numvars
[varname, xtype, varDimIDs, varAtts] = netcdf.inqVar(ncid,jj-1);
varname=varname % display name of the variables

% strcmp is a logical function used here for finding each on for the following variables: lon, lat, sst, time
% we know the names from varname and now want to pull out the values for these variables
C_lon=strcmp(varname,'longitude'); 
C_lat=strcmp(varname,'latitude'); 
C_time=strcmp(varname,'time');
C_var=strcmp(varname,'p55.162');

if C_lon == 1
% Get variable ID of the first variable, given its name.
varid = netcdf.inqVarID(ncid,varname);
% assign that variable values to x
% so x will be longitudes (in degrees East)
x = netcdf.getVar(ncid,varid);
end

if C_lat == 1
varid = netcdf.inqVarID(ncid,varname);
% y will be latitudes in degrees North
y = netcdf.getVar(ncid,varid);
end

if C_time == 1 
varid = netcdf.inqVarID(ncid,varname);
time = netcdf.getVar(ncid,varid);
end

if C_var == 1 
varid = netcdf.inqVarID(ncid,varname);
% var will have sst values
var = netcdf.getVar(ncid,varid);
end

end

% closing the nc file
netcdf.close(ncid)

tt=size(var);
% set missing values to NaN
var(var == -32767) = NaN;
% set sst to deg C (and scale it accoring to the parameters given in the file)
var=var.*scale_factor + add_offset -273.15;  

% var is 3-D matrix, here I am reshaping it into 2-D
% where rows are grid points, columns are time points
for kj=1:tt(3)
var1=squeeze(var(:,:,kj));
vard(1:tt(1)*tt(2),kj)=reshape(var1,tt(1)*tt(2),1);
end

% converting integer values to double precision 
vard=double(vard);

% lets plot the sst over the whole domain for Jan 1979
var_plot=reshape(vard(:,1),tt(1),tt(2));
load coast
figure(1); clf
h=imagesc(x,y,var_plot');
colormap jet
colorbar
axis xy
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 75]); 
xlim([190 280])
hold off

% plot VIWV for first twelve days of October 2015
figure(2); clf
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard(:,ii),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
%caxis([22 31])
if ii==12
colorbar
end
title(['October ', num2str(ii)])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% pre-processing: 
% remove seasonal cycle and apply 3-month running mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% removing 6 months in 2015 (so the timeseries is for 1979-2014)
% vard_cut=vard(:,1:end-6);
% 
% % calculate seasonal cycle for each grid point
% for i=1:length(vard(:,1))
% var_seasonal(i,:)=nanmean((reshape(vard_cut(i,:),12,36))');
% end
% 
% % copy-paste this cycle for all 36 years
% var_seasonal_all=repmat(var_seasonal,1,36);
% % add Jan-Jun for 2015 so the seasonal cycle timeseries is the same lenght as original timeseries
% var_seasonal_all=[var_seasonal_all var_seasonal(:,1:6)];
% % remove the seasonal cycle from the original series 
% vard_new=vard-var_seasonal_all;

vard_new=vard;
% apply 3-day running mean
window=3;
for i=1:length(vard(1,:))-window+1;
    runmean(i,:)=nanmean(vard_new(:,i:i+window-1)');
end
vard_runmean=runmean';


% change time to days since Oct 1 2015
time        = 1:length(time);
timeMonth   = [1/31:1/31:0.9999,... %October
               1:1/30:1.9999,... %November
               2:1/31:2.9999,... %Decomber
               3:1/31:3.9999,... %Jan
               4:1/29:4.9999,... %Feb
               5:1/31:5.9999,... %March
               6:1/30:6.9999,... %April
               7:1/31:7.1614];   %May
% remove first and last time dates because of averaging
time = time(2:end-1);   timeMonth = timeMonth(2:end-1);

% % lets plot all these pre-processing steps for one grid point (grid point number 50):
% figure;
% subplot(4,1,1)
% plot(timen,vard(50,:),'b-');
% xlabel('Time (yr)');
% ylabel('SST (^oC)');
% title('original time series for one grid point')
% xlim([1979 2015.6]);
% ylim([24 30]);
% 
% subplot(4,1,2)
% plot(timen,var_seasonal_all(50,:),'b-');
% xlabel('Time (yr)');
% ylabel('SST (^oC)');
% title('average seasonal for the grid point')
% xlim([1979 2015.6]);
% ylim([24 30])
% 
% subplot(4,1,3)
% plot(timen,vard_new(50,:),'b-');
% xlabel('Time (yr)');
% ylabel('SST anomaly (^oC)');
% title('residual (original-seasonal cycle)')
% xlim([1979 2015.6]);
% 
% subplot(4,1,4)
% plot(time,vard_runmean(50,:),'b-');
% xlabel('Time (yr)');
% ylabel('SST anomaly (^oC)');
% title('3-month running mean of the residual')
% xlim([1979 2015.6]);


% plot the 3-day running mean for all
% grid point for start of Oct
figure(3); clf
for ii=1:12
subplot(3,4,ii)
var_plot=reshape(vard_runmean(:,ii),tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
%caxis([-3 3])
if ii==12
colorbar
end
title(['runmean October ' num2str(ii)])
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PCA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% preparing the data for PCA analysis
data=vard_runmean';

% PCA:
[eigenvectors,PCs,eigenvalues]=princomp(data);

% contribution of each mode to total variance:
variance=eigenvalues./sum(eigenvalues);

figure(2); clf 
subplot(2,1,1)
plot(1:length(variance),variance,'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

subplot(2,1,2)
plot(1:20,variance(1:20),'*-');
ylabel('ratio of total variance explained by mode');
xlabel('mode number');

% plot first 4 modes (eigenvectors and PCs)
figure(3); 
for i=1:4
subplot(4,2,2*i-1)
var01=eigenvectors(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
title(['eigenvector var=',num2str(variance(i),'%2.2f')]);
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 75]); 
xlim([190 280])

subplot(4,2,2*i)
plot(timeMonth,PCs(:,i),'b-'); hold on
plot([0 max(xlim)],[0 0],'k--')
xlabel('time');
title(['PC',num2str(i)]);
ax = gca;   ax.XTick = find(timeMonth==floor(timeMonth));     ax.XTickLabel = months;
end

% plot just the first mode
figure(4); clf
for i=1
subplot(2,1,2*i-1)
var01=eigenvectors(:,i);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
axis xy
colormap jet
colorbar;
title(['eigenvector var=',num2str(variance(i),'%2.2f')]);
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 75]); 
xlim([190 280])

subplot(2,1,2*i)
plot(time,PCs(:,i),'b-')
xlabel('time');
title(['PC',num2str(i)]);
ax = gca;   ax.XTick = find(timeMonth==floor(timeMonth));     ax.XTickLabel = months;
end



%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOM algorithm 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
var01(1:tt(1)*tt(2))=NaN;
% initilizing SOM
% chose size of the map for SOM 
ny_som=4; nx_som=5;
en=ny_som*nx_som;

data=vard_runmean';
data=double(data);

msize=[ny_som nx_som];
% performing linear initialization of nodes
display('initialization')
sMap=som_lininit(data,'msize',msize,'hexa','sheet');

% training SOM
display('training')
[sM,sT] = som_batchtrain(sMap,data,'ep','hexa','sheet','radius',[3 1],'trainlen',200); 

% calulating quantization error
[q,t]=som_quality(sM,data)

% calulating hits (frequencies) of occurences of each pattern, for each seasn
hi=som_hits(sM,data);
hi=100*hi/sum(hi);

% plotting the sMap and SOM
index=[1:en];
index=reshape(index,ny_som,nx_som);
index=index';
index=reshape(index,1,nx_som*ny_som);

figure(4); clf
for i=1:en
subplot(ny_som,nx_som,i);
var01=sM.codebook(index(i),:);
var_plot=reshape(var01,tt(1),tt(2));
h=imagesc(x,y,var_plot');
colormap jet
if i == en
colorbar;
end
caxis([min(min(sM.codebook)) max(max(sM.codebook))])
axis xy  
set(h,'alphadata',~isnan(var_plot'))  
hold on
plot(long,lat,'k');
plot(long+360.,lat,'k');
ylim([45 75]); 
xlim([190 280])

% plot(xx_coast_NARR,yy_coast_NARR,'k-','LineWidth',1);
% hold on
% plot(x_dom,y_dom,'-.','Color',[0.6 0.6 0.6],'LineWidth',1);
% xlim([min(x_dom),max(x_dom)]);
% ylim([min(y_dom),max(y_dom)]); 
title([num2str(index(i)) ', f=' num2str(hi(index(i)),'%2.1f')])
axis off 
end

% find the best-matching units from the map for the given vectors
bmus=som_bmus(sM,data);
bmus_new=bmus;

% this is bmus w/o 29 Feb for leap years
bmusn=bmus_new(~isnan(bmus_new));

% plot bmus in 2D
x1=1:365;
y1=1979:2010;
figure(5); clf
imagesc(1,timeMonth,reshape(bmusn,length(bmusn),1));
ax = gca;   ax.YTickLabel = months; ax.XTickLabel = [];
colormap jet
colorbar;
%xlabel('Time (yr)');
%ylabel('Day of year');

% plot colored nodes in SOM
if en == 20
imi=[1:4; 5:8; 9:12; 13:16; 17:20]; 
elseif en == 35
imi=[1:5; 6:10; 11:15; 16:20; 21:25; 26:30; 31:35];
else
imi=[1:3; 4:6; 7:9; 10:12];
end

figure(6); clf
imagesc(imi');
colormap jet
colorbar;
axis off