% Snow density import

[snowpit_density, snowpit_text, snowpit_raw] = xlsread('summary_densitydata.xlsx','snowpit','A2:D11');

[SWEtube_density, SWEtube_text, SWEtube_raw] = xlsread('summary_densitydata.xlsx','SWEtube','A2:V34');

%Zigzag SWE locations

index = [1,3,5,8,10,13,26,28,30,32];
zigzagSWE = SWEtube_density(index',:);

for i = 1:size(zigzagSWE,1)
    for j = 1:size(zigzagSWE,2)
        if zigzagSWE(i,j)==0
           zigzagSWE(i,j-1:j)=nan; 
        elseif zigzagSWE(i,j)==1
           zigzagSWE(i,j)=nan; 
        end
    end
end

zigzagSWE = nanmean(zigzagSWE(:,1:15),2)/1000;

zigzagSWE = [SWEtube_text(index',:), num2cell(zigzagSWE)];

clear index i j 