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
    
%% Make table of coeffs
!R CMD BATCH BMS_matlab.R

load R2mat.mat

tableCol = {'Coefficient','SD', 'PIP','PSP'};
tableRow = {'aspect','elevation','northness','profileCurv','slope','tangentCurve','Sx','centreD','intercept'};

BMSG4.uniform = table(G4coeffs.Post_Mean,  G4coeffs.Post_SD, G4coeffs.PIP, G4coeffs.Cond_Pos_Sign,...
                    'RowNames', tableRow, 'VariableNames',tableCol);
BMSG4.fixed = table(G4coeffs.Post_Mean_1, G4coeffs.Post_SD_1, G4coeffs.PIP_1, G4coeffs.Cond_Pos_Sign_1,...
                    'RowNames', tableRow, 'VariableNames',tableCol);
BMSG4.random = table(G4coeffs.Post_Mean_2,G4coeffs.Post_SD_2,  G4coeffs.PIP_2, G4coeffs.Cond_Pos_Sign_2,...
                    'RowNames', tableRow, 'VariableNames',tableCol);                       
BMSG4.variable = table(G4coeffs.Post_Mean_3, G4coeffs.Post_SD_3, G4coeffs.PIP_3, G4coeffs.Cond_Pos_Sign_3,...
                    'RowNames', tableRow, 'VariableNames',tableCol);        
BMSG4.mcmc = table(G4coeffs.Post_Mean_4, G4coeffs.Post_SD_4, G4coeffs.PIP_4, G4coeffs.Cond_Pos_Sign_4,...
                    'RowNames', tableRow, 'VariableNames',tableCol);        
        
BMSG2.uniform = table(G2coeffs.Post_Mean,  G2coeffs.Post_SD, G2coeffs.PIP, G2coeffs.Cond_Pos_Sign,...
                    'RowNames', tableRow, 'VariableNames',tableCol);
BMSG2.fixed = table(G2coeffs.Post_Mean_1, G2coeffs.Post_SD_1, G2coeffs.PIP_1, G2coeffs.Cond_Pos_Sign_1,...
                    'RowNames', tableRow, 'VariableNames',tableCol);
BMSG2.random = table(G2coeffs.Post_Mean_2,G2coeffs.Post_SD_2,  G2coeffs.PIP_2, G2coeffs.Cond_Pos_Sign_2,...
                    'RowNames', tableRow, 'VariableNames',tableCol);                       
BMSG2.variable = table(G2coeffs.Post_Mean_3, G2coeffs.Post_SD_3, G2coeffs.PIP_3, G2coeffs.Cond_Pos_Sign_3,...
                    'RowNames', tableRow, 'VariableNames',tableCol);        
BMSG2.mcmc = table(G2coeffs.Post_Mean_4, G2coeffs.Post_SD_4, G2coeffs.PIP_4, G2coeffs.Cond_Pos_Sign_4,...
                    'RowNames', tableRow, 'VariableNames',tableCol); 
                
BMSG13.uniform = table(G13coeffs.Post_Mean,  G13coeffs.Post_SD, G13coeffs.PIP, G13coeffs.Cond_Pos_Sign,...
                    'RowNames', tableRow, 'VariableNames',tableCol);
BMSG13.fixed = table(G13coeffs.Post_Mean_1, G13coeffs.Post_SD_1, G13coeffs.PIP_1, G13coeffs.Cond_Pos_Sign_1,...
                    'RowNames', tableRow, 'VariableNames',tableCol);
BMSG13.random = table(G13coeffs.Post_Mean_2,G13coeffs.Post_SD_2,  G13coeffs.PIP_2, G13coeffs.Cond_Pos_Sign_2,...
                    'RowNames', tableRow, 'VariableNames',tableCol);                       
BMSG13.variable = table(G13coeffs.Post_Mean_3, G13coeffs.Post_SD_3, G13coeffs.PIP_3, G13coeffs.Cond_Pos_Sign_3,...
                    'RowNames', tableRow, 'VariableNames',tableCol);        
BMSG13.mcmc = table(G13coeffs.Post_Mean_4, G13coeffs.Post_SD_4, G13coeffs.PIP_4, G13coeffs.Cond_Pos_Sign_4,...
                    'RowNames', tableRow, 'VariableNames',tableCol);  
                
                
%% Ploting coeffs
clf
heads    = fieldnames(BMSG4);
rowNames = BMSG4.uniform.Properties.RowNames;
    for j = 1:length(heads)
        param           = char(heads(j));
        s1 = subplot(1,3,1); title('G4')
            plot(1:9,BMSG4.(param){:,1},'o','MarkerSize',10); hold on
            set(gca,'XTick', 1:length(rowNames), 'XTickLabel',rowNames); rotateticklabel(s1,45);
        s2 = subplot(1,3,2); title('G2')
            plot(1:9,BMSG2.(param){:,1},'o','MarkerSize',10); hold on
            set(gca,'XTick', 1:length(rowNames), 'XTickLabel',rowNames); rotateticklabel(s2,45);            
        s3 = subplot(1,3,3); title('G13')
            plot(1:9,BMSG13.(param){:,1},'o','MarkerSize',10); hold on
            set(gca,'XTick', 1:length(rowNames), 'XTickLabel',rowNames); rotateticklabel(s3,45);
    end
          legend(heads,'Location','best')

        fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',13)
        fig.PaperUnits = 'inches'; fig.PaperPosition = [0 0 14 7];
filename = 'BMScoeff_compare';
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')        
                
                
                
                