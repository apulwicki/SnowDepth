function [ d ] = variogramAlex(z, lag, maxlag, titletext)
%Calculates the variogam values for input data
%   This function calculates the variance between two points and their
%   distance apart. It then bins the pairs of points into specified lags
%   and calculates the average variance (semi-variance). The results are
%   returned as a strctured array with semi-variance, binned distance, and
%   number of pairs. Results are plotted. 

%   Lag is speficied by the user (lag) and the lag tolerance is set as half 
%   the lag. For example, if the lag is set at 10 then the distances will be
%   binned between 0-10, 10-20 etc, and the points plotted at 5, 15 etc.
%   Maximum lag distance can be specified but can be set to 'default' where
%   it is determined to be half the domain distance. 

%       Alexandra Pulwicki  Created: August 2016
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Calculating semi-variance

% Check formating of input data - converts to [value x y]
    if isstruct(z) %identifies transect data (fat format)
        data = [nanmean(z(5).depth(:,1:4),2), z(5).depth(:,6:7)]; %calculates mean depth
    else %identifes zigzag data
        data = z;
    end
    
% Calculate distance and variance between pairs of points
    variotest = zeros(size(data,1)*size(data,1),1); %set max size of matrix
    for i = 1:size(data,1) %over all input data
        z = data(i,1); x = data(i,2); y= data(i,3); %set the reference point
        for j = i:size(data,1) %for all new pairs with the reference point
            variotest(i*j,1) = EuclideanDistance(x,y,data(j,2),data(j,3)); %calculate euclidean distance between points
            variotest(i*j,2) = (z-data(j,1))^2; %calculate variance
        end
    end
    variotest = variotest(variotest(:,1)~=0,:); %remove zero values 
    
% Set maxlag if not specified    
    if strcmp(maxlag,'default') %when 'default' is chosen
        maxlag = ceil(max(variotest(:,1))/lag/2)*lag; %maxlag is half the maximum domain distance (rounded up to have integer number of lags)
    end
    
% Binning data based on lag
    variodata = zeros(maxlag/lag,2); %set initial matrix size
    placehold = 0:lag/2:maxlag; %contains bin ranges and labels
    for n = 2:2:length(placehold) %scan through bin ranges
        range = [placehold(1,n-1), placehold(1,n+1)]; %set range for the bin
        index = intersect(find(variotest(:,1) > range(1,1)), find(variotest(:,1) <= range(1,2))); %find where distances are within bin range
        variance = sum(variotest(index,2))/(2*length(index)); %calculate semi-variance of this set of pairs
        variodata(n/2,1:3) = [placehold(1,n), variance, length(index)]; %create matrix with the bin label, semi-variance, and number of pairs used
    end

% Create structure with output data
    d = struct('distance', variodata(:,1), 'val', variodata(:,2), 'num', variodata(:,3));

%% Plotting data

figure(1)

% Plot of variance vs lag distance (bin label)
subplot(3,1,1:2)
        h=d.distance;
        gammaexp = d.val;
        a0 = 15; % initial value: range 
        c0 = 0.1; % initial value: sill 
    [~,~,~,SS] = variogramfit(h,gammaexp,a0,c0,[],'solver','fminsearchbnd',...
                           'nugget',0,'plotit',true,...
                           'model','spherical'); hold on;   %Spherical fit
    [~,~,~,SG] = variogramfit(h,gammaexp,a0,c0,[],'solver','fminsearchbnd',...
                           'nugget',0,'plotit',true,...
                           'model','gaussian');             % Gaussian fit        
        str = {['R^2_S = ',num2str(round(SS.Rs,3))],['R^2_G = ',num2str(round(SG.Rs,3))]};
        t = annotation('textbox',[.17 .7 .2 .2],'string',str,'FitBoxToText','on');
        s = t.FontSize; t.FontSize = 14;
        title(titletext)        
        
% Coarsely binned histogram of # of pairs
subplot(3,1,3)
        %Binning of num
        lagbar = round(maxlag/lag);
        group_dist = (0:lagbar:round(max(d.distance),-1))';
        group_num = zeros(size(group_dist,1)-1,1);
        for j = 1:size(group_dist,1)-1
            index = intersect(find(d.distance>group_dist(j,1)), find(d.distance<group_dist(j+1,1)));
            group_num(j,1) = sum(d.num(index,1));
        end
    bar(mean([group_dist(1:end-1), group_dist(2:end)],2),group_num,'BarWidth', 1)
        xlabel('Lag'); ylabel('# Pairs');

% Inset plot of number of pairs        
        axes('Position',[.71 .46 .15 .13])
        box on
        plot(d.distance,d.num)
        axis([0 d.distance(end,1) 0 max(d.num)])
        ylabel('# pairs'); xlabel('lag');

end

