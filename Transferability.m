%% Transfer Coeffs between glaciers
fields = fieldnames(topo_sampled.G4);

for g = 1:3; 
    glacier = options.glacier{g}; 
for d = 1:8
    den = DenOpt{d};
    coeffs = fullLR.(den).coeff{1:7,g};
    
    for gg = 1:3; 
        glacier2 = options.glacier{gg}; 
    Transfer.(glacier).(den).(glacier2) = repmat(fullLR.(den).coeff{8,gg},options.mapsize(gg,:));
    Transfer.(glacier).(den).(glacier2)(options.mapNaN.(glacier2)) = NaN;
    
    for f = 1:length(f);
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

