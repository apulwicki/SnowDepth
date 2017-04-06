function [ rmse ] = RMSE( inputD, SWE )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global options

    for g = 1:3
       glacier = char(options.glacier(g)); 
       
       rmse.(glacier) = sqrt(sum((SWE(g).swe-inputD.(glacier)).^2)/length(SWE(g).swe));
    end

end

