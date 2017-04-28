function values = variofitAlex(d, titletext, PlotOption)
%Fits theoretical variogram curve to variogram data and plots results
%   This function fits a spherical curve to the variogram calculated by
%   variogramAlex.m. The fit is weighted by the nuber of point pairs that
%   make up each point on the variogram. It also computes the goodness of 
%   fit (Rsquared) and outputs the range, nugget and sill from the 
%   theroetical fit. Results are plotted.
%
%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Spherical Fitting
% Set initial guesses for range, nugget, and sill 
    range = median(d.meanDist); nugget = min(d.val); sill = max(d.val); 

% Set fit type to spherical
    f = fittype('nugget + ( sill*( 1.5*(h/range) - 0.5*(h/range).^3).*(h <= range) + sill*(h>range))',...
        'independent','h');

% Calculate best range, nugget and sill parameters for fit
    %Fit is weighted by number of pair points that make up each variogram
    %point
    [myfit.spherical, gof.spherical]= fit(d.binCentre,d.val,f, 'StartPoint',[nugget,sill,range],'Weights',d.num);
        nugget = round(myfit.spherical.nugget);
% Set nugget = 0 if fit nugget < 0
    if nugget < 0
        f = fittype('0 + ( sill*( 1.5*(h/range) - 0.5*(h/range).^3).*(h <= range) + (sill+nugget)*(h>range))',...
                'independent','h');
        [myfit.spherical, gof.spherical]= fit(d.binCentre,d.val,f, 'StartPoint',[sill,range],'Weights',d.num);
        nugget = 0;
    end
    
% Return value for fit parameters
    range = round(myfit.spherical.range);  sill = round(myfit.spherical.sill + nugget);
    
% Assign values for output    
    values.spherical.range = range;       values.spherical.nugget = nugget;       values.spherical.sill = sill;

%% Exponential Fitting
% Set initial guesses for range, nugget, and sill 
    range = median(d.meanDist); nugget = min(d.val); sill = max(d.val); 

% Set fit type to spherical
    f = fittype('nugget + sill*(1-exp(-h/range))', 'independent','h');

% Calculate best range, nugget and sill parameters for fit
    %Fit is weighted by number of pair points that make up each variogram
    %point
    [myfit.exponential, gof.exponential]= fit(d.binCentre,d.val,f, 'StartPoint',[nugget,sill,range],'Weights',d.num);
        nugget = round(myfit.exponential.nugget);
% Set nugget = 0 if fit nugget < 0
    if nugget < 0
        f = fittype('0 + sill*(1-exp(-h/range))', 'independent','h');
        [myfit.exponential, gof.exponential]= fit(d.binCentre,d.val,f, 'StartPoint',[sill,range],'Weights',d.num);
        nugget = 0;
    end
    
% Return value for fit parameters
    range = round(myfit.exponential.range);  sill = round(myfit.exponential.sill + nugget);
    
% Assign values for output    
    values.exponential.range = range;       values.exponential.nugget = nugget;       values.exponential.sill = sill;
    
%% Gaussian Fitting
% Set initial guesses for range, nugget, and sill 
    range = median(d.meanDist); nugget = min(d.val); sill = max(d.val); 

% Set fit type to linear
    f = fittype('nugget + ( sill*( 1 - exp(-3*h.^2/range^2)))',...
        'independent','h');

% Calculate best range, nugget and sill parameters for fit
    %Fit is weighted by number of pair points that make up each variogram
    %point
    [myfit.gaussian, gof.gaussian]= fit(d.binCentre,d.val,f, 'StartPoint',[nugget,sill,range],'Weights',d.num);
        nugget = round(myfit.gaussian.nugget);
% Set nugget = 0 if fit nugget < 0
    if nugget < 0
    f = fittype('nugget + ( sill*( 1 - exp(-3*h.^2/range^2)))',...
        'independent','h');
        [myfit.gaussian, gof.gaussian]= fit(d.binCentre,d.val,f, 'StartPoint',[sill,range],'Weights',d.num);
        nugget = 0;
    end
    
% Return value for fit parameters
    range = round(myfit.gaussian.range);  sill = round(myfit.gaussian.sill + nugget);
    
% Assign values for output    
    values.gaussian.range = range;       values.gaussian.nugget = nugget;       values.gaussian.sill = sill;
    
%% Linear Fitting
% Set initial guesses for range, nugget, and sill 
    range = median(d.meanDist); nugget = min(d.val); sill = max(d.val); 
    slope = regress(d.val,[ones(size(d.meanDist)), d.meanDist]); slope = slope(2,1);

% Set fit type to linear
    f = fittype('nugget + ( slope*h.*(h < range) + sill.*(h >= range))',...
        'independent','h');

% Calculate best range, nugget and sill parameters for fit
    %Fit is weighted by number of pair points that make up each variogram
    %point
    [myfit.linear, gof.linear]= fit(d.binCentre,d.val,f, 'StartPoint',[nugget,sill,range,slope],'Weights',d.num);
        nugget = round(myfit.linear.nugget);
% Set nugget = 0 if fit nugget < 0
    if nugget < 0
    f = fittype('nugget + ( slope*h*(h < range) + sill*(h => range))',...
        'independent','h');
        [myfit.linear, gof.linear]= fit(d.binCentre,d.val,f, 'StartPoint',[sill,range,slope],'Weights',d.num);
        nugget = 0;
    end
    
% Return value for fit parameters
    range = round(myfit.linear.range);  sill = round(myfit.linear.sill + nugget);
    
% Assign values for output    
    values.linear.range = range;       values.linear.nugget = nugget;       values.linear.sill = sill;
                
%% Plotting data

if PlotOption ==1

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
        
elseif PlotOption ==2

 %Mean fit values    
 nugget = mean([myfit.spherical.nugget, myfit.exponential.nugget, myfit.gaussian.nugget, myfit.linear.nugget]);
 range  = mean([myfit.spherical.range,  myfit.exponential.range,  myfit.gaussian.range,  myfit.linear.range]);
 sill   = mean([myfit.spherical.sill,   myfit.exponential.sill,   myfit.gaussian.sill,   myfit.linear.sill]);

 % Plot of variance vs lag distance (bin label)
    plot(d.meanDist,d.val,'o'); hold on
    plot(d.binCentre,d.val,'+'); hold on
    plot(myfit.spherical);
    plot(myfit.exponential);
    plot(myfit.linear);
    plot(myfit.gaussian);
        ylim([0 max(d.val)])
        hLeg = legend('a','b');
        set(hLeg,'visible','off')
        
        str = {['R^2_S = ',num2str(round(gof.spherical.rsquare,2))],...
               ['R^2_E = ',num2str(round(gof.exponential.rsquare,2))],...
               ['R^2_G = ',num2str(round(gof.gaussian.rsquare,2))],...
               ['R^2_L = ',num2str(round(gof.linear.rsquare,2))],...
               ['nugget = ',num2str(round(nugget,1))], ...
               ['range = ',num2str(round(range,1))],...
               ['sill = ',num2str(round(sill,1))]};
        t = annotation('textbox',[.68 .3 .2 .2],'string',str,'FitBoxToText','on');
        s = t.FontSize; t.FontSize = 12;
        title(titletext)    
        xlabel('Lag (m)'); ylabel('Semi-variance')
end
end