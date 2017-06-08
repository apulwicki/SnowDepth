function [ SWE, topo_sampled ] = ObsInCell( SWE, topo_sampled )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
global options

if length(SWE) == 3;
    for g = 1:3
    glacier = char(options.glacier(g));
    sameG   = SWE(g).cellN;

     %TOPO PARAMS
    fields = fieldnames(topo_sampled.(glacier));
        for f = 1:length(fields)
                [A1, I]  = sort(sameG);
        param = char(fields(f));
        topo_sampled.(glacier).(param)      = topo_sampled.(glacier).(param)(I,:);
        
        T = diff(A1)==0;    T1 = [0;T(1:end-1)]; T = any([T,T1],2);
        sameG_not = unique(A1(T));
        for i = length(sameG_not):-1:1
           ind = find(sameG_not(i)==A1);
           topo_sampled.(glacier).(param)(ind(1,1),:) = ...
               nanmean(topo_sampled.(glacier).(param)(ind(:,1),:));                      
           topo_sampled.(glacier).(param)(ind(2:end,1),:) = [];           
           A1(ind(2:end,1))     = [];
        end
        end
        
     %SWE
    fields = fieldnames(SWE(g));
        for f = 1:length(fields)
                [A1, I]  = sort(sameG);
        param   = char(fields(f));
        data    = SWE(g).(param)(I,:); %sort everyone
        EE      = options.E.(glacier)(I,:); %sort everyone
        NN      = options.N.(glacier)(I,:); %sort everyone
        
        T = diff(A1)==0;    T1 = [0;T(1:end-1)]; T = any([T,T1],2);
        sameG_not = unique(A1(T));
        for i = length(sameG_not):-1:1
           ind = find(sameG_not(i)==A1);
           if isnumeric(data)
           data(ind(1,1),1)      = mean(data(ind));
                if strcmp(param,'swe')
                    data(ind(1,1),2) = std(data(ind)); 
                     data(data(:,2)==0,2) = NaN;
                    data(ind,3) = (data(ind)-mean(data(ind)));%/std(data(ind));
                     data(data(:,3)==0,3) = NaN;
                    data(ind,4) = length(ind);
                     data(data(:,4)==0,4) = NaN;
                    data(ind(1,1),5) = floor(EE(ind(1,1)));
                     data(data(:,5)==0,5) = floor(EE(data(:,5)==0)); 
                    data(ind(1,1),6) = floor(NN(ind(1,1)));
                     data(data(:,6)==0,6) = floor(NN(data(:,6)==0)); 
                end
           end
           data(ind(2:end,1),:)    = [];
           A1(ind(2:end,1))     = [];
        end
        SWE(g).(param) = data;     
        end
        
    end
    
    for g = 1:3
       SWE(g).cellstd   = SWE(g).swe(:,2);
       SWE(g).standard  = SWE(g).swe(:,3);
       SWE(g).numObs    = SWE(g).swe(:,4); 
        SWE(g).numObs(isnan(SWE(g).numObs)) = 1;
       SWE(g).EN        = SWE(g).swe(:,5:6);
       SWE(g).swe(:,2:end) = [];       
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
elseif length(SWE) == 1;
    
    for g = 1:3
    glacier = char(options.glacier(g));
    sameG   = SWE.(glacier)(:,4);    
    
     %TOPO PARAMS
    fields = fieldnames(topo_sampled.(glacier));
        for f = 1:length(fields)
                [A1, I]  = sort(sameG);
        param = char(fields(f));
        topo_sampled.(glacier).(param)      = topo_sampled.(glacier).(param)(I,:);
        
        T = diff(A1)==0;    T1 = [0;T(1:end-1)]; T = any([T,T1],2);
        sameG_not = unique(A1(T));
        for i = length(sameG_not):-1:1
           ind = find(sameG_not(i)==A1);
           topo_sampled.(glacier).(param)(ind(2:end,1),:) = [];           
           A1(ind(2:end,1))     = [];
        end
        end
        
     %SWE    
    [A1, I]  = sort(sameG);
    data     = SWE.(glacier)(I,:);   %sort everyone
    T   = diff(A1)==0;    T1 = [0;T(1:end-1)]; T = any([T,T1],2);
    sameG_not = unique(A1(T));
    for i = length(sameG_not):-1:1
       ind = find(sameG_not(i)==A1);
       data(ind(1,1),1:3)    = mean(data(ind,1:3));
       data(ind(2:end,1),:)  = [];
       A1(ind(2:end,1))      = [];
    end
    SWE.(glacier) = data;
    end
   
    
end    
end

