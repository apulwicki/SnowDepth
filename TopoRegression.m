%%%%%%%%%%%%%%%%%%
% Regression between topographic parameters from SPOT5 DEM and snow depth
% measurements along transects
%%%%%%%%%%%%%%%%%%

%% Import snow depth values


%% Import Topo Params


%% Plots
elev = DEM_G4(1:3:length(DEM_G4));
slope = Slope_G4(1:3:length(Slope_G4));
aspect = Aspect_G4(1:3:length(Aspect_G4));
Sx = Sx_d100_a110_G4(1:3:length(Sx_d100_a110_G4));

mean_depth = (Depth1 + Depth2 + Depth3)/3;
mean_depth = mean_depth(1:end-1);
mean_swe = mean_depth*mean([322.3,368.3,334.7])/1000;

% fig=gcf;
% set(findall(fig,'-property','FontSize'),'FontSize',40) 

figure(1)
plot(elev, mean_swe, '.','MarkerSize',15)
xlab = xlabel('Elevation (m a.s.l.)');
ylab = ylabel('SWE (cm w.e.)');
set(xlab, 'fontsize', 30)
set(ylab, 'fontsize', 30)
l1=lsline;
set(l1,'Color',[1    0.6  0],'LineWidth',3)
fit_elev = fitlm(elev,mean_swe);
text(2300,20,strcat('R^2=',num2str(round(fit_elev.Rsquared.Ordinary,3))),'fontsize', 30)
set(gca,'fontsize', 30)
hgexport(figure(1),'G04_elev')

figure(2)
plot(slope, mean_swe, '.','MarkerSize',15)
xlab = xlabel('Slope (degrees)');
ylab = ylabel('SWE (cm w.e.)');
set(xlab, 'fontsize', 30)
set(ylab, 'fontsize', 30)
l2=lsline;
set(l2,'Color',[1    0.6  0],'LineWidth',3)
fit_slope = fitlm(slope,mean_swe);
text(60,20,strcat('R^2=',num2str(round(fit_slope.Rsquared.Ordinary,3))),'fontsize', 30)
set(gca,'fontsize', 30)
hgexport(figure(2),'G04_slope')

figure(3)
plot(aspect, mean_swe, '.','MarkerSize',15)
xlab = xlabel('Aspect (degrees)');
ylab = ylabel('SWE (cm w.e.)');
set(xlab, 'fontsize', 30)
set(ylab, 'fontsize', 30)
xlim([0 360])
l3 = lsline;
set(l3,'Color',[1    0.6  0],'LineWidth',3)
fit_aspect = fitlm(aspect,mean_swe);
text(233,20,strcat('R^2=',num2str(round(fit_aspect.Rsquared.Ordinary,3))),'fontsize', 30)
set(gca,'fontsize', 30)
hgexport(figure(3),'G04_aspect')

figure(4)
plot(Sx, mean_swe, '.','MarkerSize',15)
xlab = xlabel('Redistribution parameter');
ylab = ylabel('SWE (cm w.e.)');
set(xlab, 'fontsize', 30)
set(ylab, 'fontsize', 30)
xlim([-45 25])
l4 = lsline;
set(l4,'Color',[1    0.6  0],'LineWidth',3)
fit_Sx = fitlm(Sx,mean_swe);
text(-2,20,strcat('R^2=',num2str(round(fit_Sx.Rsquared.Ordinary,3))),'fontsize', 30)
set(gca,'fontsize', 30)
hgexport(figure(4),'G04_Sx')


% figure(2)
% plot(elev,gps_elev(1:216),'.')
% xlabel('DEM Elevation')
% ylabel('GPS Elevation')
% lsline
% fit_gpselev = fitlm(elev,gps_elev(1:216));
% %text(2050,2550,strcat('y=',num2str(fit_gpselev.Coefficients(2,1)),'x+',num2str(fit_gpselev.Coefficients.)))
% text(2050,2500,strcat('R^2=',num2str(round(fit_gpselev.Rsquared.Ordinary,2))))
