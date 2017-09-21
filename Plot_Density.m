% Plotting Snow Density Data
%       This script used imported snow density data and plots density from
%       snowpit and SWE tube. Also plots against elevation

%       Inputs:         None
%       Other scripts:  SnowDensity_Import.m
%       Outputs:        None (just saves the plots)

%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Snowpit vs SWE tube density (all glaciers)
clf
    x = cell2mat(Density.pitANDtube(:,7)); %Snowpit
    y = cell2mat(Density.pitANDtube(:,2)); %SWE tube
    errory = [cell2mat(Density.pitANDtube(:,2))-cell2mat(Density.pitANDtube(:,4)), ...
        cell2mat(Density.pitANDtube(:,5))-cell2mat(Density.pitANDtube(:,2))]; %min and max tube
    errorx = [cell2mat(Density.pitANDtube(:,7))-cell2mat(Density.pitANDtube(:,9)), ...
        cell2mat(Density.pitANDtube(:,10))-cell2mat(Density.pitANDtube(:,7))]; %min and max SP
errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','s',...
   'MarkerFaceColor','k','LineWidth',1,'MarkerSize',7); hold on
    
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y); %fit is not significant (p = 0.1)
% plot(x,yfit,'r')

%     dim = [0.15,0.79,0.11,0.11];
%     str = {strcat('y= ',num2str(round(P(1),2)),'x+ ',num2str(round(P(2)))), ...
%         strcat('R^2= ',num2str(round(LM.Rsquared.Ordinary,2)))} ;
%     annotation('textbox',dim,'String', str,'FitBoxToText','on')
%     %text(x+2, y+3, Density.pitANDtube(:,1))
    axis([220 400 220 400])
    xlabel('Snowpit-derived integrated density (kg m^{-3})')
    ylabel('Federal Sampler-derived density (kg m^{-3})')
    line = refline(1,0);
        line.Color = 'k'; line.LineStyle = '--'; hold on
     
    %Label points
    for i = 1:length(x)
        strG = char(Density.pitANDtube(i,1));
        text(x(i,1)-1, y(i,1)-4, strG,'HorizontalAlignment','right');
    end
    
    fig=gcf;    set(findall(fig,'-property','FontSize'),'FontSize',11) 
    ax = gca; ax.XTick = [220:40:400]; ax.YTick = [220:40:400];

    saveFIG('SnowpitVsSWEtube_all');

    %clear P LM x y index* i j count yfit str dim filename error* ans
    
%% SWEtube density vs elevation (all glaciers)
clf; clc;
    y = cell2mat(Density.tube(1:7,2)); %Glacier 4
    x = cell2mat(Density.tube(1:7,9));
    errory = [cell2mat(Density.tube(1:7,2))-cell2mat(Density.tube(1:7,4)),...
        cell2mat(Density.tube(1:7,5))-cell2mat(Density.tube(1:7,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit1, gof1] = fit(x,y,'poly1');
h1 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(1,:),'LineWidth',1,'MarkerSize',11); hold on
l1 = plot(SPfit1,'k'); l1.Color = options.RGB(1,:); hold on

    y = cell2mat(Density.tube(8:14,2)); %Glacier 2
    x = cell2mat(Density.tube(8:14,9));
    errory = [cell2mat(Density.tube(8:14,2))-cell2mat(Density.tube(8:14,4)),...
        cell2mat(Density.tube(8:14,5))-cell2mat(Density.tube(8:14,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit2, gof2] = fit(x,y,'poly1');
h2 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(2,:),'LineWidth',1,'MarkerSize',11); hold on
l2 = plot(SPfit2,'k'); l2.Color = options.RGB(2,:); hold on

    y = cell2mat(Density.tube(15:end,2)); %Glacier 13
    x = cell2mat(Density.tube(15:end,9));
    errory = [cell2mat(Density.tube(15:end,2))-cell2mat(Density.tube(15:end,4)),...
        cell2mat(Density.tube(15:end,5))-cell2mat(Density.tube(15:end,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit3, gof3] = fit(x,y,'poly1');
h3 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(3,:),'LineWidth',1,'MarkerSize',11); hold on
l3 = plot(SPfit3,'k'); l3.Color = options.RGB(3,:); hold on


    ylabel('SWE tube density (kg m^{-3})')
    xlabel('Elevation(m)')
    legend([h1(1) l1 h2(1) l2 h3(1) l3], {'G04',['G04 fit R^2=',num2str(round(gof1.rsquare,2))], ...
        'G02',['G02 fit R^2=',num2str(round(gof2.rsquare,2))],...
        'G13',['G13 fit R^2=',num2str(round(gof3.rsquare,2))]}, 'Location','southeast')
    fig=gcf;set(findall(fig,'-property','FontSize'),'FontSize',16)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 7.3 7];
    
 saveFIG('ElevationVsSWEtube_all');
   
    clear P LM x y index* i j count yfit str dim filename error

    
display('G04 FS'); 
display([num2str(round(SPfit1.p1,2)),'x + ',num2str(round(SPfit1.p2,0))])

display('G02 FS'); 
display([num2str(round(SPfit2.p1,2)),'x + ',num2str(round(SPfit2.p2,0))])

display('G13 FS');
display([num2str(round(SPfit3.p1,2)),'x + ',num2str(round(SPfit3.p2,0))])
    
%% Snowpit density vs elevation (all glaciers)
clf; clc

figure(3)
    y = cell2mat(Density.snowpit(1:4,2));
    x = cell2mat(Density.snowpit(1:4,5));
    errory = [cell2mat(Density.snowpit(1:4,2))-cell2mat(Density.snowpit(1:4,6)),...
        cell2mat(Density.snowpit(1:4,7))-cell2mat(Density.snowpit(1:4,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit2, gof2] = fit(x,y,'poly1');
h2 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(2,:),'LineWidth',1,'MarkerSize',11); hold on
l2 = plot(SPfit2,'k'); l2.Color = options.RGB(2,:); hold on

    y = cell2mat(Density.snowpit(5:7,2));
    x = cell2mat(Density.snowpit(5:7,5));
    errory = [cell2mat(Density.snowpit(5:7,2))-cell2mat(Density.snowpit(5:7,6)),...
        cell2mat(Density.snowpit(5:7,7))-cell2mat(Density.snowpit(5:7,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit1, gof1] = fit(x,y,'poly1');
h1 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(1,:),'LineWidth',1,'MarkerSize',11); hold on   
l1 = plot(SPfit1,'b'); l1.Color = options.RGB(1,:); hold on

    y = cell2mat(Density.snowpit(8:end,2));
    x = cell2mat(Density.snowpit(8:end,5));
    errory = [cell2mat(Density.snowpit(8:end,2))-cell2mat(Density.snowpit(8:end,6)),...
        cell2mat(Density.snowpit(8:end,7))-cell2mat(Density.snowpit(8:end,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit3, gof3] = fit(x,y,'poly1');
h3 = errorbarxy(x,y,errorx,errory,'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(3,:),'LineWidth',1,'MarkerSize',11); hold on
l3 = plot(SPfit3,'r'); l3.Color = options.RGB(3,:); hold on

    ylabel('Snowpit density (kg m^{-3})')
    xlabel('Elevation(m)')
    legend([h1(1) l1 h2(1) l2 h3(1) l3], {'G04',['G04 fit R^2=',num2str(round(gof1.rsquare,2))], ...
            'G02',['G02 fit R^2=',num2str(round(gof2.rsquare,2))],...
            'G13',['G13 fit R^2>0.99']},'Location','northeast')
    
    fig = gcf; set(findall(fig,'-property','FontSize'),'FontSize',12)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 7 6.5];
    
    saveFIG('ElevationVsSnowpit_all')

    clear P LM x y index* i j count yfit str dim filename error

    
display('G04 SP'); 
display([num2str(round(SPfit1.p1,2)),'x + ',num2str(round(SPfit1.p2,0))])
    
display('G02 SP'); 
display([num2str(round(SPfit2.p1,2)),'x + ',num2str(round(SPfit2.p2,0))])

display('G13 SP');
display([num2str(round(SPfit3.p1,2)),'x + ',num2str(round(SPfit3.p2,0))])
%% Depth vs Density - all

figure(1); clf
%SWE tube
tubeI = [1,26;27,53;54,106];
pitI = [5 7; 1 4; 8 10];
for i = 1:length(tubeI)
   t.(['p',num2str(i)]) = plot(cell2mat(Density.SWEdepth(tubeI(i,1):tubeI(i,2),3)),... %depth
       cell2mat(Density.SWEdepth(tubeI(i,1):tubeI(i,2),4)),'.',... %density
       'Color',options.RGB(i,:),'MarkerSize',20); hold on; 
   p.(['p',num2str(i)]) = plot(cell2mat(Density.snowpit(pitI(i,1):pitI(i,2),8)),... %depth
       cell2mat(Density.snowpit(pitI(i,1):pitI(i,2),2)),'o',... %density
       'Color',options.RGB(i,:),'MarkerSize',18); hold on;
end
    ylabel('Density (kg m^{-3})'); xlabel('Snow depth (cm)');
    legend([t.p1, t.p2, t.p3, p.p1, p.p2, p.p3],...
        {'Fed. Sampler - Glacier 4','Fed. Sampler - Glacier 2','Fed. Sampler - Glacier 13',...
        'Snowpit - Glacier 4','Snowpit - Glacier 2','Snowpit - Glacier 13'},'Location','southeast')
     
%fit for tube
y = cell2mat(Density.SWEdepth(:,4)); x = cell2mat(Density.SWEdepth(:,3));
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y); %fit is significant
    plot(x,yfit,'k')
        dim = [0.65,0.7,0.11,0.11];
        str = strcat('R_{FS}^2= ',num2str(round(LM.Rsquared.Ordinary,2)));
        annotation('textbox',dim,'String', str,'EdgeColor','none')
%fit for SP
y = cell2mat(Density.snowpit(:,2)); x = cell2mat(Density.snowpit(:,8));
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y); %fit is significant
    plot(x,yfit,':k')
        dim = [0.75,0.45,0.11,0.11];
        str = strcat('R_{SP}^2= ',num2str(round(LM.Rsquared.Ordinary,2)));
        annotation('textbox',dim,'String', str,'EdgeColor','none')
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',12)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 7 6.5];
        
        saveFIG('DepthDensity_SWEonly');
        
        clear coeff dim fig i LM p P pitI RegressC str t tubeI x y yfit
        
%% Depth vs Density - SP and SWE tube

figure(1); clf
%SWE tube
tubeI = [1,26;27,53;54,106];
for i = 1:length(tubeI)
   plot(cell2mat(Density.SWEdepth(tubeI(i,1):tubeI(i,2),3)),... %depth
       cell2mat(Density.SWEdepth(tubeI(i,1):tubeI(i,2),4)),'.',... %density
       'Color',options.RGB(i,:),'MarkerSize',20); hold on; 
end
    ylabel('Density (kg m^{-3})'); xlabel('Snow depth (cm)');
    legend('Glacier 4','Glacier 2','Glacier 13','Location','southeast')
%fit for tube
    y = cell2mat(Density.SWEdepth(:,4)); x = cell2mat(Density.SWEdepth(:,3));
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y); %fit is significant
    plot(x,yfit,'k')
    dim = [0.71,0.27,0.11,0.11];
    str = strcat('R^2= ',num2str(round(LM.Rsquared.Ordinary,2)));
    annotation('textbox',dim,'String', str,'FitBoxToText','on')
            fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',12)
            saveFIG('DepthDensity_tube');

    
    
figure(2); clf
%SP 
pitI = [5 7; 1 4; 8 10];
for i = 1:length(tubeI)
  plot(cell2mat(Density.snowpit(pitI(i,1):pitI(i,2),8)),... %depth
       cell2mat(Density.snowpit(pitI(i,1):pitI(i,2),2)),'.',... %density
       'Color',options.RGB(i,:),'MarkerSize',20); hold on;
end
    ylabel('Density (kg m^{-3})'); xlabel('Snow depth (cm)');    
    legend('Glacier 4','Glacier 2','Glacier 13','Location','southeast')
%fit for SP
    y = cell2mat(Density.snowpit(:,2)); x = cell2mat(Density.snowpit(:,8));
    P = polyfit(x,y,1); yfit = P(1)*x+P(2);
    LM = fitlm(x,y); %fit is significant
    plot(x,yfit,'k')
    dim = [0.71,0.27,0.11,0.11];
    str = strcat('R^2= ',num2str(round(LM.Rsquared.Ordinary,2)));
    annotation('textbox',dim,'String', str,'FitBoxToText','on')    
            fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',12)
            saveFIG('DepthDensity_SP');
           
        clear coeff dim fig i LM p P pitI RegressC str t tubeI x y yfit        
        
%% Depth vs Elevation
% Swe tube
% depth = cell2mat(Density.tube(:,10));
% elev = cell2mat(Density.tube(:,9));
% index = [1,7; 8,14; 15,31];
%     figure
% for i = 1:3
%    plot(depth(index(i,1):index(i,2)), elev(index(i,1):index(i,2)), '.',...
%        'MarkerSize', 15, 'Color', options.RGB(i,:)); hold on
% end
%     xlabel('Depth (cm)'); ylabel('Elevation (m a.s.l.)');
%     legend('Glacier 4','Glacier 2','Glacier 13')
    
% All data
% run OPTIONS.m
% options.ZZ = 2; %exclude zigzags
% run MAIN
    figure(1); clf
for j = 1:3
    depth = SWE(j).depth;
    elev = SWE(j).utm(:,3);
    i = ['n',num2str(j)];
    [f.(i), g.(i)] = fit(elev,depth,'poly1');
    subplot(1,3,j)
        h = plot(elev,depth,'.','Color',options.RGB(j,:),'MarkerSize',13); hold on
        p = plot(f.(i)); hold on
        set(p,'Color',options.RGB(j,:)); set(p, 'LineWidth',1.5);
        
        b = gca; legend(b,'off');
        dim = [b.Position(1)+0.01 b.Position(2)+.5 .3 .3];
        if j ==1;
            annotation('textbox',dim,'String', {char(options.glacier(j)),['R^2<0.01']},'FitBoxToText','on')
        else        
            annotation('textbox',dim,'String', {char(options.glacier(j)),['R^2=',num2str(g.(i).rsquare,'%.2f')]},'FitBoxToText','on')
        end
    axis([2000 2600 0 350 ])
    ylabel('Depth (cm)'); xlabel('Elevation (m a.s.l.)');
end 
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',12) 
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 13 5];

    saveFIG('DepthElevation')
%% Basic stats

% Snowpit
clc
display('G04 SP'); range = 5:7;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G02 SP'); range = 1:4;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G13 SP'); range = 8:10;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);
display('All SP'); range = 1:10;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.snowpit(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.snowpit(range,2))))),...
        ', n = ',num2str(length(range))]);

% SWE tube
display('')
display('G04 tube'); range = 1:7;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G02 tube'); range = 8:14;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
display('G13 tube'); range = 15:31;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
display('All tube'); range = 1:31;
display(['mean = ',num2str(round(nanmean(cell2mat(Density.tube(range,2))))),...
        ', std = ',num2str(round(nanstd(cell2mat(Density.tube(range,2))))),...
        ', n = ',num2str(length(range))]);
    
%% Linear regression
clc

%SP
display('G04 SP'); range = 5:7;
    SPfit = fit(cell2mat(Density.snowpit(range,5)), cell2mat(Density.snowpit(range,2)),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G02 SP'); range = 1:4;
    SPfit = fit(cell2mat(Density.snowpit(range,5)), cell2mat(Density.snowpit(range,2)),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])

display('G13 SP'); range = 8:10;
    SPfit = fit(cell2mat(Density.snowpit(range,5)), cell2mat(Density.snowpit(range,2)),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])

%SP
display('')
display('G04 tube'); range = 1:7;
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G02 tube'); range = 8:14;
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])
    
display('G13 tube'); range = [15,17,18,20:33];
    SPfit = fit(cell2mat(Density.tube(range,22)), nanmean(cell2mat(Density.tube(range,2)),2),'poly1');
    display([num2str(round(SPfit.p1,3)),'x + ',num2str(round(SPfit.p2,1))])

%% Density Range

%Snowpit
    %Name, ref, min ,max, range
disp = [Density.snowpit(:,1), num2cell(round(cell2mat(Density.snowpit(:,2)))),num2cell(round(cell2mat(Density.snowpit(:,6)))), num2cell(round(cell2mat(Density.snowpit(:,7))))];
disp = [disp, num2cell(cell2mat(disp(:,4))-cell2mat(disp(:,3)))];

%Tube
disp2 = [Density.tube(:,1), Density.tube(:,6), Density.tube(:,2), Density.tube(:,4:5)...
    num2cell(round((cell2mat(Density.tube(:,5))-cell2mat(Density.tube(:,4)))/cell2mat(Density.tube(:,2))*100))];
disp2(:,6:11) = [];



%% DETRENDED SWEtube density vs elevation (all glaciers)
options.TubeDensity     = 2;
run MAIN


clf;
    y = cell2mat(Density.tube(1:7,2)); %Glacier 4
    x = cell2mat(Density.tube(1:7,9));
    errory = [cell2mat(Density.tube(1:7,2))-cell2mat(Density.tube(1:7,4)),...
        cell2mat(Density.tube(1:7,5))-cell2mat(Density.tube(1:7,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit1, gof1] = fit(x,y,'poly1');
h1 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(1,:),'LineWidth',1,'MarkerSize',11); hold on
l1 = plot(SPfit1,'k'); l1.Color = options.RGB(1,:); hold on

    y = cell2mat(Density.tube(8:14,2)); %Glacier 2
    x = cell2mat(Density.tube(8:14,9));
    errory = [cell2mat(Density.tube(8:14,2))-cell2mat(Density.tube(8:14,4)),...
        cell2mat(Density.tube(8:14,5))-cell2mat(Density.tube(8:14,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit2, gof2] = fit(x,y,'poly1');
h2 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(2,:),'LineWidth',1,'MarkerSize',11); hold on
l2 = plot(SPfit2,'k'); l2.Color = options.RGB(2,:); hold on

    y = cell2mat(Density.tube(15:end,2)); %Glacier 13
    x = cell2mat(Density.tube(15:end,9));
    errory = [cell2mat(Density.tube(15:end,2))-cell2mat(Density.tube(15:end,4)),...
        cell2mat(Density.tube(15:end,5))-cell2mat(Density.tube(15:end,2))];
    errorx = [zeros(length(y),1),zeros(length(y),1)];
    [SPfit3, gof3] = fit(x,y,'poly1');
h3 = errorbarxy(x,y,errorx(:,2),errory(:,2),errorx(:,1),errory(:,1),'Color','k','LineStyle','none','Marker','o',...
   'MarkerFaceColor',options.RGB(3,:),'LineWidth',1,'MarkerSize',11); hold on
l3 = plot(SPfit3,'k'); l3.Color = options.RGB(3,:); hold on


    ylabel('SWE tube density (kg m^{-3})')
    xlabel('Elevation(m)')
    legend([h1(1) l1 h2(1) l2 h3(1) l3], {'G04',['G04 fit R^2=',num2str(round(gof1.rsquare,2))], ...
        'G02',['G02 fit R^2<0.01'],...
        'G13',['G13 fit R^2=',num2str(round(gof3.rsquare,2))]}, 'Location','best')
    fig=gcf;set(findall(fig,'-property','FontSize'),'FontSize',16)
    fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 7.3 7];
    
 saveFIG('ElevationVsSWEtube_all_DETRENDED');
   
    clear P LM x y index* i j count yfit str dim filename error
    
    
    %% FS-derived vs SP-derived density (fancy like in paper I)
    
        clear
load TopoSWE.mat allDepth options Density
%SP vs FS
    den     = [nan(6,1), cell2mat(Density.pitANDtube(:,2:10))];
    SP      = den(:,7);
    FS      = den(:,2);
    errorSP = [den(:,7)-den(:,9),den(:,10)-den(:,7)];   %min and max SP
    errorFS = [den(:,2)-den(:,4),den(:,5)-den(:,2)];    %min and max FS
    markerC = [options.RGB(1,:);options.RGB(1,:);...
               options.RGB(2,:);options.RGB(2,:);...
               options.RGB(3,:);options.RGB(3,:)]; 
    markerS = {'s','o','s','o','^','s'};

for i = 1:length(SP)
errorbarxy(SP(i),FS(i),errorSP(i,2),errorFS(i,2),errorSP(i,1),errorFS(i,1),...
                'Color','k','LineStyle','none','Marker',markerS{i},...
                'MarkerFaceColor',markerC(i,:),'LineWidth',1,'MarkerSize',8,...
                'MarkerEdgeColor','none'); hold on
end
    axis([220 400 220 400])
    grid on
    line = refline(1,0);
        line.Color = 'k'; line.LineStyle = '--'; hold on
    xlabel('SP-derived density (kg m^{-3})')
    ylabel('FS-derived density (kg m^{-3})')
     
    %Label points
    labels = {'G4\_USP';'G4\_LSP';'G2\_USP';'G2\_LSP';'G13\_ASP';'G13\_USP'};
    for g = 1:length(SP)
        strG = labels{g,1};
        if g ==1
        text(SP(g,1)+45, FS(g,1)-6, strG,...
            'HorizontalAlignment','right','FontSize',10,'FontName','Arial');            
        else
        text(SP(g,1)-3, FS(g,1)+7, strG,...
            'HorizontalAlignment','right','FontSize',10,'FontName','Arial');
        end
    end
    ax = gca; ax.XTick = 220:40:400; ax.YTick = 220:40:400;
        

saveFIG_IGS('SPvsFS',1,8)