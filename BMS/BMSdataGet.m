%Getting Sx data

[d300, d300text] = xlsread('d300_h0.xlsx','sampling_d300_h0','B1:BU3936');
[d200, d200text] = xlsread('d200_h0.xlsx','sampling_d200_h0','A1:BT3936');
[d100, d100text] = xlsread('d100_h0.xlsx','sampling_d100_h0','A1:BT3936');
    d300text = strcat('d300',d300text);
    d200text = strcat('d200',d200text);
    d100text = strcat('d100',d100text);
   
%     d100 = [d100text; num2cell(d100)]; 
%     d200 = [d200text; num2cell(d200)];
%     d300 = [d300text; num2cell(d300)];
    
    DDD = [d100,d200,d300]; 
     
    options.ZZ = 1; MAIN
    dataDDD = [SWE(1).pattern~='ZZ';SWE(2).pattern~='ZZ';SWE(3).pattern~='ZZ'] ;
    dataDDD = [logical(1);dataDDD];
    options.ZZ = 2; MAIN
    
    DG4     = DDD(dataDDD,:);   DG4 = DG4(1:length(SWE(1).swe),:);
    DG2     = DDD(dataDDD,:);   DG2 = DG2(length(SWE(1).swe)+1:(length(SWE(1).swe)+ length(SWE(2).swe)),:);
    DG13    = DDD(dataDDD,:);  DG13 = DG13((length(SWE(1).swe)+ length(SWE(2).swe)+2):end,:);
    
    DG4     = [SWE(1).swe,DG4];
    DG2     = [SWE(2).swe,DG2];
    DG13    = [SWE(3).swe,DG13];
    
%     DG4     = [num2cell(SWE(1).swe),DG4];
%     DG2     = [num2cell(SWE(2).swe),DG2];
%     DG13    = [num2cell(SWE(3).swe),DG13];
    
    save SxBMS.mat DG4 DG2 DG13
    
    Sx.G4 = DG4;    Sx.G2 = DG2;    Sx.G13 = DG13;  Sx.headers = ['swe',d100text, d200text, d300text];
    save SxBMS.mat Sx
    
    clear d1* d2* d3*
    
    %% Get full topo
    heads = fieldnames(topo_full.G4);
    for i = 1:3
        glacier = char(options.glacier(i));
        for j = 1:length(heads)
            param = char(heads(j));
            topoF.(glacier).(param) = topo_full.(glacier).(param)(:);
        end
    end    
    
    G4topo = topoF.G4; G2topo = topoF.G2; G13topo = topoF.G13; 
    
save fullTopo.mat G4topo G2topo  G13topo       
        
        
        
        
        
        