% Load coordinate vectors of all possible survey points

trials = 1000;     % Set number of random selections of n samples


run OPTIONS
file_path = '/Users/Alexandra/Documents/SFU/Data/SnowDepth/';
% file_path = '/home/glaciology1/Documents/QGIS/Donjek_Glaciers/Sampling/';


for g = 1:3;    glacier = options.glacier{g};
minE = min(options.rig.(glacier)(:,1));
minN = min(options.rig.(glacier)(:,2));

    nE = options.mapsize(g,2);            
    nN = options.mapsize(g,1);
    
    utmGridE.(glacier) = repmat([1:nE]*40+minE,nN,1);   
    utmGridN.(glacier) = repmat([nN:-1:1]'*40+minN,1,nE); 
end

pattern.Circle = csvread([file_path,'CellNum_Circle.csv'],1,0);
pattern.Centreline = csvread([file_path,'CellNum_Centreline.csv'],1,0);
pattern.CentreTransect = csvread([file_path,'CellNum_Transverse.csv'],1,0);
pattern.Hourglass = csvread([file_path,'CellNum_Hourglass.csv'],1,0);
pattern.HourCircle = csvread([file_path,'CellNum_HourglassCircle.csv'],1,0);

namesP = fieldnames(pattern);
for p = 5%1:length(namesP)
        [~,ia,~] = unique(pattern.(namesP{p})(:,3));
        ia = sort(ia);
    pattern.(namesP{p}) = pattern.(namesP{p})(ia,:);
    
 clear I*
for g = 1:3;    glacier = options.glacier{g};
        I_E(:,1) = pattern.(namesP{p})(:,1) > min(utmGridE.(glacier)(:));
        I_E(:,2) = pattern.(namesP{p})(:,1) < max(utmGridE.(glacier)(:));
    I_E = all(I_E,2);
        I_N(:,1) = pattern.(namesP{p})(:,2) > min(utmGridN.(glacier)(:));
        I_N(:,2) = pattern.(namesP{p})(:,2) < max(utmGridN.(glacier)(:));
    I_N = all(I_N,2);
    I = all([I_E I_N],2);
    
    E = pattern.(namesP{p})(I,1);
    N = pattern.(namesP{p})(I,2);
    
    fullUTM.(namesP{p}).(glacier) = pattern.(namesP{p})(I,1:2);
    
    m = length(N);
    
for n = 8%:30              % Set number of samples
    display(['Pattern: ',namesP{p},', Glacier: ',glacier,', n: ',num2str(n)])
    dmetric = 0;
    clear dset Nset Eset Nbest Ebest
for k=1:trials
    pp=randperm(m);
    Nset = N(pp(1:n));
    Eset = E(pp(1:n));
    
    for i=1:n
        for j=1:n
            dset(i,j) = sqrt((Nset(i)-Nset(j)).^2 + (Eset(i) - Eset(j)).^2);
        end
    end
    
     dbest = min(min(dset + 1e6*eye(n,n)));
    if dbest > dmetric
        dmetric = dbest;
        Ebest = Eset;
        Nbest = Nset;
    end
end
    utmPattern.(namesP{p})(n).(glacier)(:,1:2) = [Ebest, Nbest]; 
end
end
end
% save('PaperII_RegularSampling_SFU.mat','fullUTM','utmPattern')
%     clear g* I* ia m* n* p* utmGrid* file_path
 %%   
plot_pattern = 'Hourglass';
n = 8;

figure(1); clf
for g = 1:3;    glacier = options.glacier{g};
    Eset = utmPattern.(plot_pattern)(n).(glacier)(:,1);
    Nset = utmPattern.(plot_pattern)(n).(glacier)(:,2);

    E = fullUTM.(plot_pattern).(glacier)(:,1);
    N = fullUTM.(plot_pattern).(glacier)(:,2);
    
subplot(1,3,g)
    plot(Eset,Nset,'o'); hold on
    plot(E,N,'.')
    axis image
    axis equal
    axis([min(E) max(E) min(N) max(N)])
end
 