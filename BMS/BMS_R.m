function [ output_args ] = BMS_R( input_args )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

%% Make Cal Val data sets

y = 1:649;

 %Choose number of runs
runs = 1000;        

 %Cross validation random number matrix
[~, cal_ind] = sort(rand(runs,length(y)),2); %create matrix of random numbers
cal_ind = cal_ind(:,1:floor(length(y)*3/4)); %Choose 3/4 for the calibration component


%% Run BMS code in R
% BMS with cross validation

for i = 1:3%runs                              %for number of runs
        cal_ind_temp    = cal_ind(i,:);         %get random numbers for choosing obs
        val_ind         = setdiff(1:length(y),cal_ind_temp); %get random numbers for validating

        sweG4   = SWE(1).swe(cal_ind_temp);
        sweG2   = SWE(2).swe(cal_ind_temp);
        sweG13  = SWE(3).swe(cal_ind_temp);

        heads   = fieldnames(topo_sampled_ns.G4);
        Vtopo.G4 = [];  Vtopo.G2 = [];  Vtopo.G13 = []; 
        for j = 1:length(heads)
            param             = char(heads(j));
            topoG4.(param)  = topo_sampled_ns.G4.(param)(cal_ind_temp);
            topoG2.(param)  = topo_sampled_ns.G2.(param)(cal_ind_temp); 
            topoG13.(param) = topo_sampled_ns.G13.(param)(cal_ind_temp);
            Vtopo.G4        = [Vtopo.G4, topo_sampled_ns.G4.(param)(val_ind)];
            Vtopo.G2        = [Vtopo.G2, topo_sampled_ns.G2.(param)(val_ind)]; 
            Vtopo.G13       = [Vtopo.G13, topo_sampled_ns.G13.(param)(val_ind)];
        end

        save mat2R.mat sweG4 sweG2 sweG13 topoG4 topoG2 topoG13
        % Run BMS code in R
        !R CMD BATCH BMS_matlab.R

        % Make table of coeffs
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

            models = fieldnames(BMSG4);

            bms.uniform(i,:)    = {BMSG4.uniform, BMSG2.uniform, BMSG13.uniform};
            bms.fixed(i,:)      = {BMSG4.fixed, BMSG2.fixed, BMSG13.fixed};
            bms.random(i,:)     = {BMSG4.random, BMSG2.random, BMSG13.random};
            bms.variable(i,:)   = {BMSG4.variable, BMSG2.variable, BMSG13.variable};
            bms.mcmc(i,:)       = {BMSG4.mcmc, BMSG2.mcmc, BMSG13.mcmc};

            
            % validation
            for g = 1:3
                glacier = char(options.glacier(g));
            for m = 1:length(models)
                mod = char(models(m));
                
                y_regress = sum(Vtopo.(glacier).*repmat(bms.uniform{i,g}{1:end-1,1}',...
                                    length(Vtopo.(glacier)),1),2) + bms.uniform{i,g}{end,1};      
                rmse.(mod).(glacier)(i,1) = sqrt(sum((SWE(g).swe(val_ind)-...
                                            y_regress).^2)/numel(y_regress));
            end
            end            
            
end

    %min rmse coefficients
    for g = 1:3
        glacier = char(options.glacier(g));
        BMSbest.(glacier) = table();
    for m = 1:length(models)
        mod = char(models(m));
        
        rmse_min    = rmse.(mod).(glacier)==min(rmse.(mod).(glacier));
        bestCoefs   = [bms.(mod){rmse_min,g}(:,1); ...
            table(min(rmse.(mod).(glacier)),'VariableNames',{'Coefficient'},'RowNames',{'rmse'})]; 
        bestCoefs.Properties.VariableNames = {mod};
        
        BMSbest.(glacier) = [BMSbest.(glacier), bestCoefs];
    end
    end    

%

%     
end

