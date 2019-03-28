function [swe_out, topo_out] = deterministic_locations(subsetSWE, subsetTOPO, n)

%% ************************************************************************
% This maximin sort-of design is deterministic and therefore unique for a 
% given n. It begins with point #1 in the list and adds subsequent points 
% based on maximizing the minimum distance between any two points. This 
% means that measurement locations are consistent for increasing n (ie the 
% sampling scheme doesn't entirely change as each new sampling location is 
% added). Various metrics were tried to create uniform sampling, including 
% max total or mean distances between all points and minimum standard 
% deviation between distances. Maximin worked the best empirically.
    
    E = subsetSWE(:,2);
    N = subsetSWE(:,3);
    SWE = subsetSWE(:,1);

q = length(N);
param = fieldnames(subsetTOPO);

% Do the first two calculations outside of the main loop 
nnn = 50;
Nbest(1) = N(nnn);
Ebest(1) = E(nnn);
SWEbest(1) = SWE(nnn);
for t = 1:length(param); TOPObest.(param{t})(1) = subsetTOPO.(param{t})(nnn);  end

Nleft = N([1:nnn-1,nnn+1:end]);
Eleft = E([1:nnn-1,nnn+1:end]);
SWEleft = SWE([1:nnn-1,nnn+1:end]);
for t = 1:length(param); TOPOleft.(param{t}) = subsetTOPO.(param{t})([1:nnn-1,nnn+1:end]);  end

clear dvec
for i=1:length(Nleft)
    dvec(i) = sqrt((Nleft(i)-Nbest(1)).^2 + (Eleft(i) - Ebest(1)).^2);
end

[furthest] = find (dvec == max(dvec));
Nbest(2) = Nleft(furthest(1));
Ebest(2) = Eleft(furthest(1));
SWEbest(2) = SWEleft(furthest(1));
for t = 1:length(param); TOPObest.(param{t})(2) = TOPOleft.(param{t})(furthest(1));  end

Nleft = Nleft([1:furthest-1,furthest+1:end]);
Eleft = Eleft([1:furthest-1,furthest+1:end]);
SWEleft = SWEleft([1:furthest-1,furthest+1:end]);
for t = 1:length(param); TOPOleft.(param{t}) = TOPOleft.(param{t})([1:furthest-1,furthest+1:end]);  end

% loop over the remainder of the samples n 
    for k=3:n  

  clear dvec
  
  for i=1:length(Nleft)  % loop over the number of unused samples
    
    for j = 1:k-1   % calculate the distance matrix to all used samples
      dvec(i,j) = sqrt((Nleft(i)-Nbest(j)).^2 + (Eleft(i) - Ebest(j)).^2);
    end
 
  end
  
  mindvec = min(dvec');
  [furthest] = find(mindvec == max(mindvec));
  
  Nbest(k) = Nleft(furthest(1));
  Ebest(k) = Eleft(furthest(1));
  SWEbest(k) = SWEleft(furthest(1));
  for t = 1:length(param); TOPObest.(param{t})(k) = TOPOleft.(param{t})(furthest(1));  end


  Nleft = Nleft([1:furthest-1,furthest+1:end]);
  Eleft = Eleft([1:furthest-1,furthest+1:end]);
  SWEleft = SWEleft([1:furthest-1,furthest+1:end]);
  for t = 1:length(param); TOPOleft.(param{t}) = TOPOleft.(param{t})([1:furthest-1,furthest+1:end]);  end

    end
    swe_out = [SWEbest', Ebest', Nbest'];
    topo_out = TOPObest;
end
