function [ SWEdata, TOPOdata ] = SortNSelect( SWEdata, TOPOdata, n )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
global options

for g = 1:3
    glacier = char(options.glacier(g));
    
     %Sort data based on cell#
    data                = SWEdata.(glacier);
    [ ~, I ]            = sort(data(:,4));
    SWEdata.(glacier)   = data(I,:);
     ff = fieldnames(TOPOdata.(glacier));
     for i = 1:length(ff);    fname = ff{i};
    TOPOdata.(glacier).(fname)  = TOPOdata.(glacier).(fname)(I,:);  end
    
     %Select number of points (evenly spaced)
    nI = floor(linspace(1,length(SWEdata.(glacier)),n));
    SWEdata.(glacier)   = data(nI,:);
         for i = 1:length(ff);    fname = ff{i};
    TOPOdata.(glacier).(fname)  = TOPOdata.(glacier).(fname)(nI,:); end
end

end

