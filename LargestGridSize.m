
%% Combine observations into larger grid sizes
S = [2:17].^2;

figure(2);clf

for s = 1:length(S)

SWEoutput = GridUpsizing( S(s) ); 

%      topoParam.G4  = NaN(options.mapsize(1,:));    topoParam.G2  = NaN(options.mapsize(2,:));    topoParam.G13 = NaN(options.mapsize(3,:));
% PlotTopoParameter(topoParam,'none', 'SWE (m w.e.)',SWEoutput,'colour', 'NOmassB')

%!!!Figure out where these werid numbers are coming from!!!
for g = 1:3; glacier = options.glacier{g};
SWEoutput.(glacier)(SWEoutput.(glacier)(:,4)>1,4) = NaN;
meanSTD.(glacier)(s) = nanmean(SWEoutput.(glacier)(:,4));
end
end

for g = 1:3; glacier = options.glacier{g};
subplot(1,3,g)
    plot(S, meanSTD.(glacier)); hold on
    %plot([0 max(S)], [ZZstd(g),ZZstd(g)],'k--')
end
%% 
    load Full.mat fullLR options
    
for g = 1:3; glacier = options.glacier{g};
num = fullLR.S2.(glacier);
ZZstd = [0.027 0.035 0.040];

S = [2:17].^2;

for s = 1:length(S)
    clear D
sizeG = S(s);
    C = sqrt(sizeG);

for i = 1:C:size(num,1)-C
   for j =  1:C:size(num,2)-C
       mini = num(i:i+C-1,j:j+C-1);
       D(i,j) = nanstd(mini(:)-nanmean(mini(:)));
   end
end
D(D==0) = [];   D(isnan(D)) =[];
GmeanSTD(s) = mean(D);

figure(1);
subplot(sqrt(length(S)), sqrt(length(S)),s)
    histogram(D)
    title(['Grid size = ',num2str(S(s))])

end

figure(2);
subplot(1,3,g)
    plot(S, GmeanSTD); hold on
    plot([0 max(S)], [ZZstd(g),ZZstd(g)],'k--')
end