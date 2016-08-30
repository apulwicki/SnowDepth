glacier = 'G04'; %select data from chosen glacier
pattern = 'all';
    %z = pulldata(data, book, glacier, person, pattern, quality, format)
z = pulldata(SD,'all',glacier,'all','all',1,'fat'); %transect data  
z1 = pulldata(SD,'Extra',glacier,'Extra','Extra',1,'fat'); %ExtraSD data from nontransect measurements

x = [z(5).depth(:,6);z1(5).depth(:,42)]; x1 = nanmax(x)-x; %convert easting to distance in m
y = [z(5).depth(:,7);z1(5).depth(:,43)]; y1 = nanmax(y)-y; %convert easting to distance in m
z = [nanmean(z(5).depth(:,1:4),2);nanmean(z1(5).depth(:,1:40),2)];

d = variogramAlex([z x1 y1], 15, 'default', glacier);

%%
close all

a = 100; n=800; sig = var(z);

f = fittype('n + ( sig*( 1.5*(h/a) - 0.5*(h/a).^3).*(h <= a) + sig*(h>a))',...
    'independent','h');

[myfit gof]= fit(d.binCentre,d.val,f, 'StartPoint',[a,n,sig],'Weights',d.num);

range = round(myfit.a); nugget = round(myfit.n); sill = round(myfit.sig);


%% Plotting data

titletext = 'test';

figure(1)

% Plot of variance vs lag distance (bin label)
subplot(3,1,1:2)
    plot(d.meanDist,d.val,'o'); hold on
    plot(d.binCentre,d.val,'+'); hold on
    plot(myfit);
        ylim([0 max(data(:,1))])
        hLeg = legend('a','b');
        set(hLeg,'visible','off')
        
        str = {['R^2_S = ',num2str(round(gof.rsquare,3))],...
                ['nugget = ',num2str(nugget)], ['range = ',num2str(range)],['sill = ',num2str(sill)]};
        t = annotation('textbox',[.25 .37 .2 .2],'string',str,'FitBoxToText','on');
        s = t.FontSize; t.FontSize = 12;
        title(titletext)        
        
% Coarsely binned histogram of # of pairs
subplot(3,1,3)
    barwidth = round(length(d.binCentre)/10);
    bins = d.binCentre(1:barwidth:end,1)-d.binCentre(1,1);
    numMtx = zeros(length(bins)-1,1);
    for i = 1:length(bins)-1
        logiEl = d.binCentre > bins(i) & d.binCentre <= bins(i+1);
        numMtx(i) = sum(d.num(logiEl));
    end   
    bar(numMtx,'BarWidth', 1)
        xlabel('Lag'); ylabel('# Pairs');

% Inset plot of number of pairs        
        axes('Position',[.71 .46 .15 .13])
        box on
        plot(d.binCentre,d.num)
        axis([0 d.binCentre(end,1) 0 max(d.num)])
        ylabel('# pairs'); xlabel('lag');

