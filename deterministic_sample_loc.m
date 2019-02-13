%% ************************************************************************
% This maximin sort-of design is deterministic and therefore unique for a 
% given n. It begins with point #1 in the list and adds subsequent points 
% based on maximizing the minimum distance between any two points. This 
% means that measurement locations are consistent for increasing n (ie the 
% sampling scheme doesn't entirely change as each new sampling location is 
% added). Various metrics were tried to create uniform sampling, including 
% max total or mean distances between all points and minimum standard 
% deviation between distances. Maximin worked the best empirically.

run OPTIONS
load('PaperII_RegularSampling.mat')

namesP = fieldnames(fullUTM);


for p = 1:length(namesP);   pattern = namesP{p};
for g=1:3;                  glacier = options.glacier{g};
    
E = fullUTM.(pattern).(glacier)(:,1);
N = fullUTM.(pattern).(glacier)(:,2);

for n = 6:30            % Set number of samples

close all
clear  i j k dset dmin dsum dvec dvecsum dbest dabsmin 
clear Ebest Nbest Eleft Nleft kkbest p Eset Nset
clear stdvec furthest

q = length(N);

% Do the first two calculations outside of the main loop 
Nbest(1) = N(1);
Ebest(1) = E(1);

Nleft = N(2:end);
Eleft = E(2:end);

clear dvec
for i=1:length(Nleft)
    dvec(i) = sqrt((Nleft(i)-Nbest(1)).^2 + (Eleft(i) - Ebest(1)).^2);
end

[furthest] = find (dvec == max(dvec));
Nbest(2) = Nleft(furthest(1));
Ebest(2) = Eleft(furthest(1));

Nleft = Nleft([1:furthest-1,furthest+1:end]);
Eleft = Eleft([1:furthest-1,furthest+1:end]);


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

  Nleft = Nleft([1:furthest-1,furthest+1:end]);
  Eleft = Eleft([1:furthest-1,furthest+1:end]);
end
    samplingSubset.(pattern)(n).(glacier)(:,1:2) = [Ebest', Nbest'];

end 
end
end

% figure
% plot(E,N,'.')
% hold on
% for i=1:k
% plot(Ebest(i),Nbest(i),'ko')
% pause
% end