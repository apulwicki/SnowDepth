function [ WBoutput ] = WBnoise( WBinput, HIGHLOW )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global options

if strcmp(HIGHLOW,'high')

    for g = 1:3;    glacier = options.glacier{g};
        WBoutput.(glacier) = WBinput.(glacier) + ...
                normrnd( 0, options.zzstd(g)*5, length(WBinput.(glacier)),1);
    end

elseif strcmp(HIGHLOW,'low')

    for g = 1:3;    glacier = options.glacier{g};
        WBoutput.(glacier) = WBinput.(glacier) + ...
                normrnd( 0, options.zzstd(g), length(WBinput.(glacier)),1);
    end    
    
end

