
%% Importing DEMs
files = dir('/home/glaciology1/Documents/Data/GlacierTopos/*.asc');
%files = dir('/home/glaciology1/Documents/Data/GlacierTopos/*.asc');
v = cell(0,0);
for i = 1:length(files)
    if i >=10 && i<=12
        filename = ['/home/glaciology1/Documents/Data/GlacierTopos/', files(i).name];
        startRow = 8;
        formatSpec = '%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%27f%f%[^\n\r]';
        fileID = fopen(filename,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', '', 'WhiteSpace', '', 'EmptyValue' ,NaN,'HeaderLines' ,startRow-1, 'ReturnOnError', false);
        fclose(fileID);
        A.data = [dataArray{1:end-1}]; A.data(A.data<-3e+38) = NaN;
        clearvars filename startRow formatSpec fileID dataArray ans;  
    elseif i >=7 && i<=9
        A = importdata(['/home/glaciology1/Documents/Data/GlacierTopos/', files(i).name],' ',6);
        A.data(A.data==0) = NaN; 
    else
        A = importdata(['/home/glaciology1/Documents/Data/GlacierTopos/', files(i).name],' ',6);
        A.data(A.data==-9999) = NaN; 
    end
    A.data(~any(~isnan(A.data), 2),:)=[];
    A.data(:,~any(~isnan(A.data), 1))=[];
    v(i,1) = cellstr(files(i).name(1:end-4)); eval([v{i,1} '= A.data;']);
end
    
for i = 1:length(v)/3
    param = char(v(i*3)); param = param(1:end-4);
    eval(['topo_full.G4.(param) = ',v{3*i-1,1},';']);
    eval(['topo_full.G2.(param) = ',v{3*i-2,1}],';');
    eval(['topo_full.G13.(param) = ',v{3*i,1}],';');
end
    clear aspect* elev* north* profil* slope* Sx* tangent* files A i param v
    
%% Import Topo Params

%Import topos
[topo, ~, ~] = xlsread('SPOT_TopoParams_noZZ.xlsx','Sheet1','A1:G2353');

%Remove zigzag values
run OPTIONS.m; options.ZZ = 2;
run MAIN.m

div = [1, length(SWE(1).swe); length(SWE(1).swe)+1, length(SWE(1).swe)+length(SWE(2).swe);...
        length(SWE(1).swe)+length(SWE(2).swe)+1, length(SWE(1).swe)+length(SWE(2).swe)+length(SWE(3).swe)];
glacier = [4,2,13];

for i = 1:3
name = ['G', num2str(glacier(i))];
    aspect           = topo(div(i,1):div(i,2),2);
    northness        = topo(div(i,1):div(i,2),3);
    profileCurve     = topo(div(i,1):div(i,2),4);
    tangentCurve     = topo(div(i,1):div(i,2),5);
    slope            = topo(div(i,1):div(i,2),6);
    elevation        = topo(div(i,1):div(i,2),7);
    
    topo_sampled.(name) = struct('aspect',aspect, 'elevation',elevation,...
        'northness',northness, 'profileCurve',profileCurve, 'slope',slope,...
        'tangentCurve',tangentCurve);
end

% Sx import & stepwise MLR
[d300, d300text] = xlsread('d300_h0.xlsx','sampling_d300_h0','B1:BU3936');
[d200, d200text] = xlsread('d200_h0.xlsx','sampling_d200_h0','A1:BT3936');
[d100, d100text] = xlsread('d100_h0.xlsx','sampling_d100_h0','A1:BT3936');
    d300text = strcat('d300',d300text);
    d200text = strcat('d200',d200text);
    d100text = strcat('d100',d100text);

for i = 1:3
    y       = SWE(i).swe;
    name    = ['G', num2str(glacier(i))];
    X       = [d100(div(i,1):div(i,2),:), d200(div(i,1):div(i,2),:), d300(div(i,1):div(i,2),:)];
    text    = [d100text, d200text, d300text];
   
    [~, ~, ~, inmodel] = stepwisefit(X,y);
    text    = text(1,inmodel);
    X1      = [ones(length(X),1), X(:,inmodel)];
    mlr_Sx     = regress(y,X1); mlr_sort = real(sort(complex(mlr_Sx)));
    
    best            = find(mlr_Sx == mlr_sort(end-1,1));
    winddir.(name)  = text(best);
    Sx              = X1(:,best);
    topo_sampled.(name).Sx = Sx;
end
        clear best d300* d200* d100* i inmodel mlr* name text X* y

        
%Distance from centreline import
run CentrelineDistance.m

%Keep a copy of non-standardized variables for range of params sampled
%plots
topo_sampled_ns = topo_sampled;

%Standardizing variables
params = fieldnames(topo_sampled.G4);
for i = 1:3
name= char(glacier(i));
    for t = 1:length(params)
    field = char(params(t));
    
    topo_sampled.(name).(field) = (topo_sampled.(name).(field)-...
        mean(topo_sampled.(name).(field)))/std(topo_sampled.(name).(field));
    end
end

    clear i name topo field j params div glacier G t X1 Y1 min_dist fields distance centreline corner
    clear aspect* elev* north* profil* slope* Sx* tangent*    

%Put zigzags back in
run OPTIONS.m; options.ZZ = 1;
run MAIN.m   
