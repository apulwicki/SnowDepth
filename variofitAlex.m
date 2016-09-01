function myfit = variofitAlex(d, titletext)

%% Spherical Fitting
range = 100; nugget = min(d.val); sill = max(d.val);

f = fittype('nugget + ( sill*( 1.5*(h/range) - 0.5*(h/range).^3).*(h <= range) + sill*(h>range))',...
    'independent','h');

[myfit, gof]= fit(d.binCentre,d.val,f, 'StartPoint',[range,nugget,sill],'Weights',d.num);

range = round(myfit.range); nugget = round(myfit.nugget); sill = round(myfit.sill+myfit.nugget);


%% Plotting data

figure(1)

% Plot of variance vs lag distance (bin label)
subplot(3,1,1:2)
    plot(d.meanDist,d.val,'o'); hold on
    plot(d.binCentre,d.val,'+'); hold on
    plot(myfit);
        ylim([0 max(d.val)])
        hLeg = legend('a','b');
        set(hLeg,'visible','off')
        
        str = {['R^2_S = ',num2str(round(gof.rsquare,3))],...
                ['nugget = ',num2str(nugget)], ['range = ',num2str(range)],['sill = ',num2str(sill)]};
        t = annotation('textbox',[.17 .7 .2 .2],'string',str,'FitBoxToText','on');
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
    bar(bins(2:end), numMtx,'BarWidth', 1)
        xlabel('Lag'); ylabel('# Pairs');

% Inset plot of number of pairs        
        axes('Position',[.71 .46 .15 .13])
        box on
        plot(d.binCentre,d.num)
        axis([0 d.binCentre(end,1) 0 max(d.num)])
        ylabel('# pairs'); xlabel('lag');

end