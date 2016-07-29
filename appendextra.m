function [ SDmean, gname ] = appendextra( SD, glacier )
%APPENDEXTRA Summary of this function goes here
%   Detailed explanation goes here
gname = glacier;

data = pulldata(SD,'all',glacier,'all','all',1,'fat');
data2 = pulldata(SD,'ZZ',glacier,'ZZ','ZZ',1,'fat');

SDmean = [nanmean(data(5).depth(:,1:4),2), nanstd(data(5).depth(:,1:4),1,2), data(5).depth(:,6), data(5).depth(:,7); ...
            nanmean(data2(5).depth(:,1:41),2), nanstd(data2(5).depth(:,1:41),1,2), data2(5).depth(:,42), data2(5).depth(:,43)];

end

