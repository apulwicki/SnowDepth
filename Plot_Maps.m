
load TopoBMS_MLR

%% Maps of topographic params for each glacier

header  = fieldnames(topo_full_ns.G4);

for r = 1%:length(header)
    param = char(header(r));
    topoParam.G4  = topo_full_ns.G4.(param);
    topoParam.G2  = topo_full_ns.G2.(param);
    topoParam.G13 = topo_full_ns.G13.(param);

    PlotTopoParameter(topoParam,param, options.topoVarsUnits(r), SWE, 'black', 'nomassB')
    
%     %Save figure
    saveFIG(['Map_',param])
end 
    clear r topoParam param header
%% SWE at sampling locations

param = 'empty';
for opt = 2%:9
topoParam.G4  = NaN(size(topo_full_ns.G4.elevation));
topoParam.G2  = NaN(size(topo_full_ns.G2.elevation));
topoParam.G13 = NaN(size(topo_full_ns.G13.elevation));
topoParam.rig = options.rig;

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', sweOPT(opt), 'colour')

    saveFIG(['SWEmap_opt',num2str(opt)])
end
%% Modelled and observed SWE

modelled = sweMLR;
type     = 'MLR';

% modelled = sweBMS;
% type     = 'BMS';

% opt      = 8;
for opt = 2:9
    topoParam.G4  = modelled(opt).G4;
    topoParam.G2  = modelled(opt).G2;
    topoParam.G13 = modelled(opt).G13;

    PlotTopoParameter(topoParam, 'modelledSWE', 'SWE (m w.e.)', SWE, 'sweONswe', 'massB')
         saveFIG([type,'map_Modelled_Observed',num2str(opt-1)])
end
    clear filename modelled opt type fig glacier g 
%% Modelled SWE Difference as %
% modelled = sweMLR;
% type = 'MLR';
modelled = sweBMS;
type = 'BMS';

    %Differencing modelled SWE
for g = 1:3
    glacier = char(options.glacier(g));
    
    for i = 2:9
    stackSWE.(glacier)(:,:,i-1)   = modelled(i).(glacier);
    end

minSWE.(glacier)  = nanmin(stackSWE.(glacier),[],3);
maxSWE.(glacier)  = nanmax(stackSWE.(glacier),[],3);
    hereNan = isnan(sweMLR(2).(glacier));
    minSWE.(glacier)(hereNan) = NaN;    maxSWE.(glacier)(hereNan) = NaN;   meanSWE.(glacier)(hereNan) = NaN; 

diffSWE.(glacier) = maxSWE.(glacier)-minSWE.(glacier);
    diffSWE.(glacier) = [nan(2,size(diffSWE.(glacier),2));diffSWE.(glacier);nan(2,size(diffSWE.(glacier),2))];
    diffSWE.(glacier) = [nan(size(diffSWE.(glacier),1),2),diffSWE.(glacier),nan(size(diffSWE.(glacier),1),2)];
diffSWE_p.(glacier) = (maxSWE.(glacier)-minSWE.(glacier))./minSWE.(glacier)*100;
    diffSWE_p.(glacier)(minSWE.(glacier)==0) = 0;    
    diffSWE_p.(glacier)(diffSWE_p.(glacier)>70) = 70;
    diffSWE_p.(glacier) = [nan(2,size(diffSWE_p.(glacier),2));diffSWE_p.(glacier);nan(2,size(diffSWE_p.(glacier),2))];
    diffSWE_p.(glacier) = [nan(size(diffSWE_p.(glacier),1),2),diffSWE_p.(glacier),nan(size(diffSWE_p.(glacier),1),2)];

end
diffSWE.rig = rig;   diffSWE_p.rig = rig; 

PlotTopoParameter(diffSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
   saveFIG([type,'_SWEdifferenceMap'])
PlotTopoParameter(diffSWE_p, 'modelledSWE', {'Predicted SWE range as ','percent of maximum SWE (%)'}, SWE, 'black')
    saveFIG([type,'_SWEdifferenceMap_percent'])

display(char(type))
for g = 1:3
    glacier = char(options.glacier(g));
    display([glacier,' ', num2str(round(nanmean(diffSWE_p.(glacier)(:))))])
end
    
    
    clear g* i maxSWE minSWE filename modelled  type diffSWE*

%% Min/Max/Mean of modelled SWE and total difference

    %Differencing modelled SWE
for g = 1:3
    glacier = char(options.glacier(g));
    
    for i = 2:9
    stackSWE.(glacier)(:,:,i-1)       = sweMLR(i).(glacier);
    stackSWE.(glacier)(:,:,2*(i-1))   = sweBMS(i).(glacier);
    end

minSWE.(glacier)    = min(stackSWE.(glacier),[],3);  
maxSWE.(glacier)    = max(stackSWE.(glacier),[],3);
meanSWE.(glacier)   = nanmean(stackSWE.(glacier),3);
    hereNan = isnan(sweMLR(2).(glacier));
    minSWE.(glacier)(hereNan) = NaN;    maxSWE.(glacier)(hereNan) = NaN;   meanSWE.(glacier)(hereNan) = NaN; 

diffSWE.(glacier)   = maxSWE.(glacier)-minSWE.(glacier);    
diffSWE_p.(glacier) = (maxSWE.(glacier)-minSWE.(glacier))./meanSWE.(glacier)*100;
    weirdsmall = diffSWE_p.(glacier)<120;
    diffSWE_p.(glacier)(weirdsmall) = 120;
    weirdsmall = diffSWE_p.(glacier)>180;
    diffSWE_p.(glacier)(weirdsmall) = 180;
end
minSWE.rig = rig;   maxSWE.rig = rig;   meanSWE.rig = rig;
sweMEAN.rig = rig;  sweMIN.rig = rig;   sweMAX.rig = rig;   sweRANGE.rig = rig;
diffSWE.rig = rig;  diffSWE_p.rig = rig;


%Plots from cell values
PlotTopoParameter(minSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black'); 
    saveFIG('SWEmin');
PlotTopoParameter(maxSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    saveFIG('SWEmax');    
PlotTopoParameter(diffSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    saveFIG('SWEdiff');
PlotTopoParameter(diffSWE_p, 'modelledSWE', 'Percent of Mean SWE (%)', SWE, 'black')
    saveFIG('SWEdiff_percentMean');

%Plots from coefficients
PlotTopoParameter(sweMEAN, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black');
    saveFIG('SWEmeanModelled');
PlotTopoParameter(sweMIN, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    saveFIG('SWEminModelled');
PlotTopoParameter(sweMAX, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    saveFIG('SWEmaxModelled');    
PlotTopoParameter(sweRANGE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    saveFIG('SWErangeModelled');    
      

clear g* i maxSWE minSWE filename modelled  type diffSWE

%% Residuals at sampling locations

param = 'empty';
for g = 1:3
    glacier = char(options.glacier(g));
topoParam.(glacier)  = NaN(size(topo_full_ns.(glacier).elevation));
topoParam.rig = rig;

resZ(g).swe = res.(glacier); 
resZ(g).utm = SWE(g).utm;
end

figure(3)
PlotTopoParameter(topoParam,param, 'BMS Residuals (m w.e.)', resZ, 'colour')
    C = cbrewer('div', 'RdYlBu', 20, 'PCHIP');
    colormap(flipud(C))
        saveFIG(['residualsMap_',method])


