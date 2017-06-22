clear alex taylor

taylor(:,1) = [571731.48;577258.98;580978.1;587346.4;591126.5;597353.2;601796.1;608101];
taylor(:,2) = [6737517.35;6733918.68;6730286.9;6730436.4;6724959.2;6730694.1;6734532.2;6736574.4];
taylor(:,3) = [2620;2640;2380;2225;2070;1915;1765;1615];
taylor(:,4) = [3.2;3.7;2.6;2.7;2.4;1.78;1.65;0.91];
taylor(:,5) = [0.407;0.409;0.394;0.364;0.388;0.385;0.390;0.342];
taylor(:,6) = [1.302;1.513;1.024;0.983;0.932;0.685;0.643;0.311];

alex(:,1) = [566453.4;570077.4;595349.9;601424.8;605031];
alex(:,2) = [6727621.2;6732429.7;6741065.2;6753607.6;6762717.2];
alex(:,3) = [2610;2730;2321;2472;2434];
alex(:,6) = [1.30;1.59;0.62;0.50;0.32];
  %  alex = alex(3:end,:);

%% Distance from Divide Regression
% Dt = sqrt((taylor(:,1)-taylor(1,1)).^2+(taylor(:,2)-taylor(1,2)).^2)/1000;
% Da = sqrt((alex(:,1)-taylor(1,1)).^2+(alex(:,2)-taylor(1,2)).^2)/1000;
Dt = sqrt((taylor(:,1)-alex(1,1)).^2+(taylor(:,2)-alex(1,2)).^2)/1000;
Da = sqrt((alex(:,1)-alex(1,1)).^2+(alex(:,2)-alex(1,2)).^2)/1000;

LMt = fitlm(Dt,taylor(:,6));
LMa = fitlm(Da,alex(:,6));

[Ft, Gt] = fit(Dt,taylor(:,6),'poly1');
    CIt = predint(Ft,Dt);
[Fa, Ga] = fit(Da,alex(:,6),'poly1');
    CIa = predint(Fa,Da);

    [Dat, Iat] = sort([Da; Dt]); AT = [alex(:,6); taylor(:,6)]; AT = AT(Iat);    
[Fat, Gat] = fit(Dat,AT,'poly1');
    CIat = predint(Fat,Dat);    
    


clf
plot(Dt,taylor(:,6),'.', 'MarkerSize',13, 'Color',rgb('DarkCyan')); hold on
    pt = plot(Ft);    set(pt,'Color',rgb('DarkCyan')); set(pt, 'LineWidth',1.5)
    plot(Dt, CIt(:,1),'--', 'Color',rgb('DarkCyan'), 'LineWidth',.2)
    plot(Dt, CIt(:,2),'--','Color',rgb('DarkCyan'), 'LineWidth',.2)

plot(Da,alex(:,6),'.', 'MarkerSize',13, 'Color',rgb('OrangeRed'))
    pa = plot(Fa);   set(pa,'Color',rgb('OrangeRed')); set(pa, 'LineWidth',1.5)
    plot(Da, CIa(:,1), '--','Color',rgb('OrangeRed'), 'LineWidth',.2)
    plot(Da, CIa(:,2),'--','Color',rgb('OrangeRed'), 'LineWidth',.2)

legend([pt,pa],'Taylor1969','Alex2016')
    xlabel('Distance from Sian Camp (km)'); ylabel('SWE (m w.e.)')
    yl = ylim;
    ylim([0, yl(2)]);
    
figure(2); clf
    pt = plot(Fat);    set(pt,'Color','k'); set(pt, 'LineWidth',1.5); hold on
%     plot(Dat, CIat(:,1),'--', 'Color','k', 'LineWidth',.2)
%     plot(Dat, CIat(:,2),'--','Color','k', 'LineWidth',.2)
L(2) = plot(Da,alex(:,6),'.', 'MarkerSize',13, 'Color',rgb('OrangeRed'));
       plot(Da,alex(:,6),'.', 'MarkerSize',13, 'Color',rgb('OrangeRed'));
L(1) = plot(Dt,taylor(:,6),'.', 'MarkerSize',13, 'Color',rgb('DarkCyan')); hold on
    xlim([min(Dat), max(Dat)]); ylim([0 2])
    grid on
        xlabel('Distance from mountain divide (km)'); ylabel('SWE (m w.e.)')
        legend(L,'Taylor (1969)','Pulwicki and Flowers (2017)', 'Location','northoutside')
saveFIG_IGS('AccumGrad',1,8.6)

%% Elevation Regression

[~, I] =sort(alex(:,3)); 
Ea = [alex(I,3), alex(I,6)];

LMt = fitlm(taylor(:,3),taylor(:,6));
LMa = fitlm(Ea(:,1),Ea(:,2));

[Ft, Gt] = fit(taylor(:,3),taylor(:,6),'poly1');
    CIt = predint(Ft,taylor(:,3));
[Fa, Ga] = fit(Ea(:,1),Ea(:,2),'poly1');
    CIa = predint(Fa,Ea(:,1));


clf
plot(taylor(:,3),taylor(:,6),'.', 'MarkerSize',13, 'Color',rgb('DarkCyan')); hold on
    pt = plot(Ft);    set(pt,'Color',rgb('DarkCyan')); set(pt, 'LineWidth',1.5)
    plot(taylor(:,3), CIt(:,1),'--', 'Color',rgb('DarkCyan'), 'LineWidth',.2)
    plot(taylor(:,3), CIt(:,2),'--','Color',rgb('DarkCyan'), 'LineWidth',.2)

plot(Ea(:,1),Ea(:,2),'.', 'MarkerSize',13, 'Color',rgb('OrangeRed'))
    pa = plot(Fa);   set(pa,'Color',rgb('OrangeRed')); set(pa, 'LineWidth',1.5)
    plot(Ea(:,1), CIa(:,1), '--','Color',rgb('OrangeRed'), 'LineWidth',.2)
    plot(Ea(:,1), CIa(:,2),'--','Color',rgb('OrangeRed'), 'LineWidth',.2)

legend([pt,pa],'Taylor1969','Alex2016')
    xlabel('Elevation (m a.s.l.)'); ylabel('SWE (m w.e.)')
    yl = ylim;
    ylim([0, yl(2)]);