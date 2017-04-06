function [ sampledP ] = SampledCell( InputGrid )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
global options
    for g = 1:3;
        glacier = char(options.glacier(g));

        data = InputGrid.(glacier);
        sampledtemp = data(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
        sampledP.(glacier) = diag(sampledtemp);
    end

end

