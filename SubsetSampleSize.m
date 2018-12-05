function [ WBout, TOPOout, UTMout ] = SubsetSampleSize( WB, TOPO, UTM, sampleSize )
%SubsetSampleSize get an evenly distributed subset of the WB and TOPO data
%of the size provided
global options

    namesP = fieldnames(WB); 

for g = 1:3;    glacier = options.glacier{g};
    for p = 1:length(namesP)

    
    %Index - evenly spaced along full data set
    I = unique(round(linspace(1,length(WB.(namesP{p}).(glacier)), sampleSize)));
    
    %WB select data
    WBout.(namesP{p}).(glacier) = WB.(namesP{p}).(glacier)(I);
    UTMout.(namesP{p}).(glacier) = UTM.(namesP{p}).(glacier)(I,:);
        
    %TOPO select data
     param = fields(TOPO.(namesP{1}).G4);
    for f = 1:length(param);    field = param{f};
    TOPOout.(namesP{p}).(glacier).(field) = TOPO.(namesP{p}).(glacier).(field)(I);
    end
    
    end
%         UTMout.Random.(glacier)(:,1) = [];
end
end

