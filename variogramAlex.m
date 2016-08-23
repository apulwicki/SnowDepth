function [ d ] = variogramAlex(z, lag, maxlag, titletext)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

if isstruct(z) %transect data
    data = [nanmean(z(5).depth(:,1:4),2), z(5).depth(:,6:7)];
else %zigzag data
    data = z;
end    
    
% Calculate distance between points (lag)

%find the variance for each pair of points
variotest = zeros(size(data,1)*size(data,1),1);
for i = 1:size(data,1)
    z = data(i,1); x = data(i,2); y= data(i,3);
    for j = 1:size(data,1)
        variotest(i*j,1) = EuclideanDistance(x,y,data(j,2),data(j,3)); 
        variotest(i*j,2) = (z-data(j,1))^2;
    end
end

if ~exist('maxlag','var')
    maxlag = round(max(variotest(:,1)),-1);
end

variodata = zeros(maxlag/lag,2);
placehold = 0:lag/2:maxlag;
for n = 3:2:length(placehold)
    range = [placehold(1,n-2), placehold(1,n)];
    index = intersect(find(variotest(:,1) > range(1,1)), find(variotest(:,1) <= range(1,2)));
    variance = sum(variotest(index,2))/(2*length(index));
    variodata((n-1)/2,1:3) = [placehold(1,n-1), variance, length(index)];
end

d = struct('distance', variodata(:,1), 'val', variodata(:,2), 'num', variodata(:,3));

%variogram fit
    figure(1)
     subplot(3,1,1:2)
        h=d.distance;
        gammaexp = d.val;
        a0 = 15; % initial value: range 
        c0 = 0.1; % initial value: sill 
        [~,~,~,SS] = variogramfit(h,gammaexp,a0,c0,[],...
                               'solver','fminsearchbnd',...
                               'nugget',0,'plotit',true,...
                               'model','spherical');
       hold on 
       [~,~,~,SG] = variogramfit(h,gammaexp,a0,c0,[],...
                               'solver','fminsearchbnd',...
                               'nugget',0,'plotit',true,...
                               'model','gaussian');                   
        str = {['R^2_S = ',num2str(round(SS.Rs,3))],['R^2_G = ',num2str(round(SG.Rs,3))]};
        t = annotation('textbox',[.17 .7 .2 .2],'string',str,'FitBoxToText','on');
        s = t.FontSize;
        t.FontSize = 14;
        title(titletext)        
        
        subplot(3,1,3)
        lagbar = lag*3;
        group_dist = (0:lagbar:round(max(d.distance),-1))';
        group_num = zeros(size(group_dist,1)-1,1);
        for j = 1:size(group_dist,1)-1
            index = intersect(find(d.distance>group_dist(j,1)), find(d.distance<group_dist(j+1,1)));
            group_num(j,1) = sum(d.num(index,1));
        end
        bar(mean([group_dist(1:end-1), group_dist(2:end)],2),group_num,'BarWidth', 1)
        xlabel('Lag'); ylabel('# Pairs');
        
        axes('Position',[.71 .46 .15 .13])
        box on
        plot(d.distance,d.num)
        axis([0 d.distance(end,1) 0 max(d.num)])
        ylabel('# pairs'); xlabel('lag');

end

