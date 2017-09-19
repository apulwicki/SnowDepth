function [ WBoutput ] = WBnoise( WBinput )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global options

    for g = 1:3;    glacier = options.glacier{g};
        WBoutput.(glacier) = WBinput.(glacier) + ...
                normrnd( 0, options.zzstd(g), length(WBinput.(glacier)),1);
    end

end

