function [ WBoutput ] = WBnoise( WBinputT, HIGHLOW )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global options

if strcmp(HIGHLOW,'high')

    for g = 1:3;    glacier = options.glacier{g};
        WBoutput.(glacier) = WBinputT.(glacier)(:,1) + ...
                normrnd( 0, options.zzstd(g)*5, size(WBinputT.(glacier),1),1);
    end

elseif strcmp(HIGHLOW,'low')

    for g = 1:3;    glacier = options.glacier{g};
        WBoutput.(glacier) = WBinputT.(glacier)(:,1) + ...
                normrnd( 0, options.zzstd(g), size(WBinputT.(glacier),1),1);
    end    
    
end

