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
    [myfit.spherical, gof.spherical]= fit(d.binCentre,d.val,f, ...
                                        'StartPoint',[nugget,sill,range],...
                                        'Weights',d.num,...
                                        'Lower',[0,0,0]);%,...
                                        %'Upper',[max(d.val),max(d.val)*5, max(d.meanDist)*5]);
              
% Assign values for output    
    values.spherical.range  = round(myfit.spherical.range);
    values.spherical.nugget = round(myfit.spherical.nugget);
    values.spherical.sill   = round(myfit.spherical.sill + myfit.spherical.nugget);
    values.spherical.gof    = gof.spherical;

%% Exponential Fitting
% Set initial guesses for range, nugget, and sill 
    range = median(d.meanDist); nugget = min(d.val); sill = max(d.val); 

% Set fit type to spherical
    f = fittype('nugget + sill*(1-exp(-h/range))', 'independent','h');

% Calculate best range, nugget and sill parameters for fit
    %Fit is weighted by number of pair points that make up each variogram
    %point
    [myfit.exponential, gof.exponential]= fit(d.binCentre,d.val,f,...
                                        'StartPoint',[nugget,sill,range],...
                                        'Weights',d.num,...
                                        'Lower',[0,0,0]);%,...
                                        %'Upper',[max(d.val),max(d.val), max(d.meanDist)]);
             
% Assign values for output    
    values.exponential.range  = round(myfit.exponential.range);
    values.exponential.nugget = round(myfit.exponential.nugget);
    values.exponential.sill   = round(myfit.exponential.sill + myfit.exponential.nugget);
    values.exponential.gof    = gof.exponential;
    
%% Gaussian Fitting
% Set initial guesses for range, nugget, and sill 
    range = median(d.meanDist); nugget = min(d.val); sill = max(d.val); 

% Set fit type to linear
    f = fittype('nugget + ( sill*( 1 - exp(-3*h.^2/range^2)))',...
        'independent','h');

% Calculate best range, nugget and sill parameters for fit
    %Fit is weighted by number of pair points that make up each variogram
    %point
    [myfit.gaussian, gof.gaussian]= fit(d.binCentre,d.val,f,...
                                        'StartPoint',[nugget,sill,range],...
                                        'Weights',d.num,...
                                        'Lower',[0,0,0]);%,...
                                        %'Upper',[max(d.val),max(d.val), max(d.meanDist)]);

% Assign values for output    
    values.gaussian.range  = round(myfit.gaussian.range);
    values.gaussian.nugget = round(myfit.gaussian.nugget);
    values.gaussian.sill   = round(myfit.gaussian.sill + myfit.gaussian.nugget);
    values.gaussian.gof    = gof.gaussian;

%% Linear Fitting
% % Set initial guesses for range, nugget, and sill 
%     range = median(d.meanDist); nugget = min(d.val); sill = max(d.val); 
%     slope = regress(d.val(1:round(end/2)),[ones(size(d.val(1:round(end/2)))), d.meanDist(1:round(end/2))]); 
%         slope = slope(2,1);
% 
% % Set fit type to linear
%     f = fittype('nugget + slope*h.*(h < range) + sill*(h >= range)',...
%         'independent','h');
% 
% % Calculate best range, nugget and sill parameters for fit
%     %Fit is weighted by number of pair points that make up each variogram
%     %point
%     [myfit.linear, gof.linear]= fit(d.binCentre,d.val,f, 'StartPoint',[nugget,sill,range,slope],'Weights',d.num);
%         nugget = round(myfit.linear.nugget);
% % Set nugget = 0 if fit nugget < 0
%     if nugget < 0
%     f = fittype('0 + ( slope*h.*(h < range) + sill*(h >= range))',...
%         'independent','h');
%         [myfit.linear, gof.linear]= fit(d.binCentre,d.val,f, 'StartPoint',[sill,range,slope],'Weights',d.num);
%         nugget = 0;
%     end
%     
% % Return value for fit parameters
%     range = round(myfit.linear.range);  sill = round(myfit.linear.sill + nugget);
%     
% % Assign values for output    
%     values.linear.range = range;       values.linear.nugget = nugget;       values.linear.sill = sill;
%                 
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
 nugget = mean([values.spherical.nugget, values.exponential.nugget, values.gaussian.nugget]);%, values.linear.nugget]);
 range  = mean([values.spherical.range,  values.exponential.range,  values.gaussian.range]);%,  values.linear.range]);
 sill   = mean([values.spherical.sill,   values.exponential.sill,   values.gaussian.sill]);%,   values.linear.sill]);

 [colormap] = cbrewer('qual', 'Set1', 3);
 
 % Plot of variance vs lag distance (bin label)
    %plot(d.meanDist,d.val,'o'); hold on
    %plot(d.binCentre,d.val,'+'); hold on
    p1 = plot(myfit.spherical); hold on;
        p1.Color = colormap(1,:); 
    p2 = plot(myfit.exponential);
         p2.Color = colormap(2,:);
    %plot(myfit.linear);
    p3 = plot(myfit.gaussian);
         p3.Color = colormap(3,:);   
        %ylim([0 max(d.val)])
        xlim([0 max(d.binCentre)])
        hLeg = legend('a','b');
        set(hLeg,'visible','off')
        
        str = {['R^2_S = ',num2str(round(gof.spherical.rsquare,2))],...
               ['R^2_E = ',num2str(round(gof.exponential.rsquare,2))],...
               ['R^2_G = ',num2str(round(gof.gaussian.rsquare,2))],...
               %['R^2_L = ',num2str(round(gof.linear.rsquare,2))],...
               ['nugget = ',num2str(round(nugget,1))], ...
               ['range = ',num2str(round(range,1))],...
               ['sill = ',num2str(round(sill,1))]};
        %t = annotation('textbox',[.68 .3 .2 .2],'string',str,'FitBoxToText','on');
        %s = t.FontSize; t.FontSize = 12;
        title(titletext)    
        xlabel('Lag (m)'); ylabel('Semi-variance')
        
end
end