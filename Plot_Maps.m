
load TopoBMS_MLR

%% Maps of topographic params for each glacier

header  = fieldnames(topo_full_ns.G4);

for r = 1:length(header)
    param = char(header(r));
    topoParam.G4  = topo_full_ns.G4.(param);
    topoParam.G2  = topo_full_ns.G2.(param);
    topoParam.G13 = topo_full_ns.G13.(param);
    topoParam.rig = rig;

    PlotTopoParameter(topoParam,param, options.topoVarsUnits(r), SWE, 'black')
    
%     %Save figure
    filename = ['Map_',param];
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
end 
    clear r topoParam param header
%% SWE at sampling locations

param = 'empty';
for opt = 2:9
    run OPTIONS
    options.DensitySWE = opt;
    run MAIN

topoParam.G4  = NaN(size(topo_full_ns.G4.elevation));
topoParam.G2  = NaN(size(topo_full_ns.G2.elevation));
topoParam.G13 = NaN(size(topo_full_ns.G13.elevation));
topoParam.rig = rig;

PlotTopoParameter(topoParam,param, 'SWE (m w.e.)', SWE, 'colour')

    filename = ['SWEmap_opt',num2str(opt)];
print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
end
%% Modelled and observed SWE

modelled = sweMLR;
type     = 'MLR';

% modelled = sweBMS;
% type     = 'BMS';

%opt      = 8;
for opt = 2:9
    topoParam.G4  = modelled(opt).G4;
    topoParam.G2  = modelled(opt).G2;
    topoParam.G13 = modelled(opt).G13;
    topoParam.rig = rig;

    PlotTopoParameter(topoParam, 'modelledSWE', 'SWE (m w.e.)', SWE, 'colour')

    %Integrated specific winter balance
       for g = 1:3; glacier = char(options.glacier(g)); 
        ISWbalance(opt).(glacier) = round(nanmean(modelled(opt).(glacier)(:)),2);
       end
       
    annotation('textbox',[.17 .22 .1 .1],'String',[num2str(ISWbalance(opt).G4, '%.2f'),' m w.e.'],'EdgeColor','none')    
    annotation('textbox',[.34 .55 .1 .1],'String',[num2str(ISWbalance(opt).G2, '%.2f'),' m w.e.'],'EdgeColor','none')    
    annotation('textbox',[.75 .53 .1 .1],'String',[num2str(ISWbalance(opt).G13, '%.2f'),' m w.e.'],'EdgeColor','none')    
    
    fig=gcf; set(findall(fig,'-property','FontSize'),'FontSize',18)
    filename = [type,'map_Modelled_Observed_Opt',num2str(opt-1)];
        print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
end
    clear filename modelled opt type
%% Modelled SWE Difference as %
modelled = sweMLR;
type = 'MLR';
% modelled = sweBMS;
% type = 'BMS';

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
diffSWE_p.(glacier) = (maxSWE.(glacier)-minSWE.(glacier))./maxSWE.(glacier)*100;

weirdsmall = diffSWE_p.(glacier)<10;
diffSWE_p.(glacier)(weirdsmall) = 10;

weirdsmall = diffSWE_p.(glacier)>100;
diffSWE_p.(glacier)(weirdsmall) = 100;

diffSWE_p.(glacier)(all(stackSWE.(glacier)==0,3)) = 0;

end

diffSWE.rig = rig;    diffSWE_p.rig = rig;


% PlotTopoParameter(diffSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
%     filename = [type,'_SWEdifferenceMap'];
%     print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
PlotTopoParameter(diffSWE_p, 'modelledSWE', {'Predicted SWE range as ','percent of maximum SWE (%)'}, SWE, 'black')
    filename = [type,'_SWEdifferenceMap_percent'];
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
     
    
    
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

%Plots from cell values
PlotTopoParameter(minSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    filename = 'SWEmin';
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
PlotTopoParameter(maxSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    filename = 'SWEmax';    
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')  
PlotTopoParameter(diffSWE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    filename = 'SWEdiff';
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')   
PlotTopoParameter(diffSWE_p, 'modelledSWE', 'Percent of Mean SWE (%)', SWE, 'black')
    filename = 'SWEdiff_percentMean';
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')     

%Plots from coefficients
PlotTopoParameter(sweMEAN, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    filename = 'SWEmeanModelled';    
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')  
PlotTopoParameter(sweMIN, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    filename = 'SWEminModelled';
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
PlotTopoParameter(sweMAX, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    filename = 'SWEmaxModelled';    
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
PlotTopoParameter(sweRANGE, 'modelledSWE', 'SWE (m w.e.)', SWE, 'black')
    filename = 'SWErangeModelled';    
    print([options.path1, filename],'-dpng','-r0'); print([options.path2, filename],'-dpng','-r0')
      

clear g* i maxSWE minSWE filename modelled  type diffSWE



