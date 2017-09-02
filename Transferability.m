load Full.mat
load TopoSWE.mat

%% Transfer Coeffs between glaciers
fields = fieldnames(topo_sampled.G4);

for g = 1:3; 
    glacier = options.glacier{g}; 
for d = 1:8
    den = options.DenOpt{d};
    coeffs = fullLR.(den).coeff{1:7,g};
    
    for gg = 1:3; 
        glacier2 = options.glacier{gg}; 
    Transfer.(glacier).(den).(glacier2) = repmat(fullLR.(den).coeff{8,g},options.mapsize(gg,:));
    Transfer.(glacier).(den).(glacier2)(options.mapNaN.(glacier2)) = NaN;
    
    for f = 1:length(fields);
       param = fields{f};
       Transfer.(glacier).(den).(glacier2) = topo_full.(glacier2).(param)*coeffs(f)...
                                        + Transfer.(glacier).(den).(glacier2);
    end
 stackTransfer.(glacier).(glacier2)(:,:,d) =  Transfer.(glacier).(den).(glacier2); 
    end
end

 for gg = 1:3; glacier2 = options.glacier{gg};
     Transfer.(glacier).mean.(glacier2) = nanmean(stackTransfer.(glacier).(glacier2),3);
 end
end

    %clear stackTransfer g*

%% Use all data
for d = 1:8;    den = options.DenOpt{d};

% Input Data
[ tempSWE, tempTOPO ] = ObsInCell(fullSWE.(den).input, topo_sampled); 

comboSWE.(den).G4 = [tempSWE.G4; tempSWE.G2; tempSWE.G13];
    comboSWE.(den).G2 = comboSWE.(den).G4; comboSWE.(den).G13 = comboSWE.(den).G4;

FF = fieldnames(tempTOPO.G4);
for f = 1:length(FF)
   param = FF{f};
   comboTOPO.G4.(param) = [tempTOPO.G4.(param); tempTOPO.G2.(param); tempTOPO.G13.(param)]; 
        comboTOPO.G2 = comboTOPO.G4; comboTOPO.G13 = comboTOPO.G4;
end

% Linear Regression 
comboLR.(den) =  LinearRegression( comboSWE.(den), comboTOPO, topo_full );

% Choose best set of coeffs
    [~,I] = min(comboLR.(den).coeff{10,:}); 
    comboLR.(den).coeff = comboLR.(den).coeff(:,I);
    coeffs = fullLR.(den).coeff{1:7,I};

% Apply coeffs
    for gg = 1:3; 
        glacier2 = options.glacier{gg}; 
    Transfer.combo.(den).(glacier2) = repmat(comboLR.(den).coeff{8,1},options.mapsize(gg,:));
    Transfer.combo.(den).(glacier2)(options.mapNaN.(glacier2)) = NaN;
    
    for f = 1:length(f);
       param = FF{f};
       Transfer.combo.(den).(glacier2) = topo_full.(glacier2).(param)*coeffs(f)...
                                        + Transfer.combo.(den).(glacier2);
    end
        Transfer.combo.(den).(glacier2)(Transfer.combo.(den).(glacier2)<0) = 0;
    end
end

% Stack and get mean of all density options
for d = 1:8
 den = options.DenOpt{d};
 stackTransfer.combo.(glacier2)(:,:,d) =  Transfer.combo.(den).(glacier2); 
end
 for gg = 1:3; glacier2 = options.glacier{gg};
     Transfer.combo.mean.(glacier2) = nanmean(stackTransfer.combo.(glacier2),3);
 end

%% PLOT -> maps of mean


figure(1); PlotTopoParameter( Transfer.G4.mean,'empty', 'SWE (m w.e.)', SWE, 'none', 'massB')
     title('Glacier 4 LR Coefficients')
     saveFIG('MapTransferabilityG4Coeffs')

figure(2); PlotTopoParameter( Transfer.G2.mean,'empty', 'SWE (m w.e.)', SWE, 'none', 'massB')
     title('Glacier 2 LR Coefficients')
     saveFIG('MapTransferabilityG2Coeffs')

figure(3); PlotTopoParameter( Transfer.G13.mean,'empty', 'SWE (m w.e.)', SWE, 'none', 'massB')
     title('Glacier 13 LR Coefficients')
     saveFIG('MapTransferabilityG13Coeffs')

figure(4); PlotTopoParameter( Transfer.combo.mean,'empty', 'SWE (m w.e.)', SWE, 'none', 'massB')
     title('All Data LR Coefficients')
     saveFIG('MapTransferabilityComboCoeffs')   
     
%% RMSE of data from each coeff set

%SWE = ObsInCell(SWE, topo_sampled);

    for d = 1:8
    den = options.DenOpt{d};
    
    G4RES(d)    = SampledCell( Transfer.G4.(den) );
    G2RES(d)    = SampledCell( Transfer.G2.(den) );
    G13RES(d)   = SampledCell( Transfer.G13.(den) );
    comboRES(d) = SampledCell( Transfer.combo.(den) );
  
    for g = 1:3
    glacier = options.glacier{g};
    
    G4RMSE(g,d)     = sqrt(mean((SWE(g).swe-G4RES(d).(glacier)).^2));
    G2RMSE(g,d)     = sqrt(mean((SWE(g).swe-G2RES(d).(glacier)).^2));
    G13RMSE(g,d)    = sqrt(mean((SWE(g).swe-G13RES(d).(glacier)).^2));
    comboRMSE(g,d)  = sqrt(mean((SWE(g).swe-comboRES(d).(glacier)).^2));
    end
    end

% Plot it
figure(1); clf

    ymax = max(max([G4RMSE;G2RMSE;G13RMSE;comboRMSE]));
C = [options.RGB; [0 67 133]/255];

data = [mean(G4RMSE,2), mean(G2RMSE,2), mean(G13RMSE,2), mean(comboRMSE,2)];
B = bar(data);
    for i = 1:4
        B(i).FaceColor = C(i,:);
    end
    legend({'G4 Coefficients', 'G2 Coefficients', 'G13 Coefficients', 'All data Coefficients'})
    ylabel('Mean RMSE (m w.e.)');
    set(gca, 'XTick', [1 2 3]);    set(gca, 'XTickLabel', {'Glacier 4' 'Glacier 2' 'Glacier 13'})   
%saveFIG('TransferabilityRMSE',13)    
    
%% Mean elev coeff and intercept 

for d = 1:8
   den = options.DenOpt{d};
   elevcoeff(d,:) = [fullLR.(den).coeff{1,:}, comboLR.(den).coeff{1,:}];
   intercept(d,:) = [fullLR.(den).coeff{8,:},comboLR.(den).coeff{8,:}];
end

elevcoeff = mean(elevcoeff);
intercept = mean(intercept);

mean(G4RMSE(:))
mean(G2RMSE(:))
mean(G13RMSE(:))
mean(comboRMSE(:))

