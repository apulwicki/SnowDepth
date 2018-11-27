% %% Coefficients for all data 
% 
% % Selecting Data from Pattern
% 
%     den = 'S2';
% 
% load TopoSWE.mat
% run OPTIONS
% clear DataObs_RMSE  
% 
% [ SWE, topo_full, STD ]    = ObsInCell(fullSWE.(den).input, topo_full);
% 
%     % Remove dc, aspect and Northness
% for g = 1:3;    glacier = options.glacier{g};
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
% end
% 
% % BMA
% swe_input = SWE;
% topo_input = topo_full;
% 
% cd BMS
% [BMSinit, BMSres] = BMS_R(swe_input, topo_input);
% cd ..
% 
% for g = 1:3;        glacier = char(options.glacier(g));
% BMS.(glacier) = BMSinit.(glacier);   
% BMS.(glacier).Properties.VariableNames = {'BMSCoefficient','BMSsemiR2','BMSunivarR2'};
% residualsBMS.(glacier) = BMSres.(glacier);
% end

% BMS_alldata = BMS
% save('PaperII_AblationArea.mat','BMS_alldata', '-append')
%% Coefficients for each pattern

% Selecting Data from Pattern
% 
%     den = 'S2';
% 
% load TopoSWE.mat
% run OPTIONS
% 
% for t = [1,3,4,5,100] %[6,1,3,4,5,100]
% if     t == 6; type = 'Circle';           subset = 'pattern';       
% elseif t == 1; type = 'Centreline';       subset = 'pattern';     
% elseif t == 3; type = 'CentreTransect';   subset = 'pattern';   
% elseif t == 4; type = 'Hourglass';        subset = 'pattern';  
% elseif t == 5; type = 'HourCircle';       subset = 'pattern';
% elseif t == 100; type = 'RandomSafe';     subset = 'random';
% end
% 
% 
% input.SWE = fullSWE.(den); input.topo_sampled = topo_sampled; 
% input.topo_sampled_ns = topo_sampled_ns;
% 
% [ subsetSWE_temp, TOPOdata_temp ] = DataSubset( subset, t, input );
% 
% [ subsetSWE_temp, TOPOdata_temp ] = ObsInCell( subsetSWE_temp, TOPOdata_temp ); 
% 
% maxN = min([length(subsetSWE_temp.G4) length(subsetSWE_temp.G2) length(subsetSWE_temp.G13)]);
% 
%     % Remove dc, aspect and Northness
% for g = 1:3;    glacier = options.glacier{g};
%     TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'centreD');
%     TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'aspect');
%     TOPOdata_temp.(glacier) = rmfield(TOPOdata_temp.(glacier),'northness');
% end
% 
% 
% 
%  for n = 5:maxN
%      display([type, ' n=',num2str(n)])
% 
% for g = 1:3;    glacier = options.glacier{g};
%     nI = randperm(maxN, n);
%     %nI = floor(linspace(1,maxN,n));
%     WBinput(n).(type).(glacier)   = subsetSWE_temp.(glacier)(nI,:);
%         ff = fieldnames(TOPOdata_temp.(glacier));
%     for i = 1:length(ff);    fname = ff{i};
%     TOPOinput(n).(type).(glacier).(fname)  = TOPOdata_temp.(glacier).(fname)(nI,:); end
% 
% end
% 
% % Linear regression (BMA)
%     swe_input = WBinput(n).(type);
%     topo_input = TOPOinput(n).(type);
%     
%     cd BMS
%     [BMSinit, BMSres] = BMS_R(swe_input, topo_input);
%     cd ..
%     
%     for gg = 1:3;        glacier = char(options.glacier(gg));
%     BMS(n).(type).(glacier) = BMSinit.(glacier);   
%     BMS(n).(type).(glacier).Properties.VariableNames = {['BMSCoefficient_', num2str(t)],...
%                                                  ['BMSsemiR2_', num2str(t)],...
%                                                  ['BMSunivarR2_', num2str(t)]};
%     residualsBMS(n).(type).(glacier) = BMSres.(glacier);    
%     end
% 
%  end
% end
% % 


%% Predict
% load PaperII_AblationArea.mat
% load TopoSWE.mat 
% run OPTIONS.m
% 
%     % Remove dc, aspect and Northness
% for g = 1:3;    glacier = options.glacier{g};
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'centreD');
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'aspect');
%     topo_full.(glacier) = rmfield(topo_full.(glacier),'northness');
% end
% 
% 
% for g = 1:3;    glacier = options.glacier{g};
%     sweBMS.(glacier) = repmat(BMS_alldata.(glacier){5,1}, options.mapsize(g,:));
%     mlrCoeff = BMS_alldata.(glacier){1:4,1};    topoCoeff = fieldnames(topo_full.G4);
%         %multiply coeffs and add them
%     for m = 1:length(mlrCoeff)
%         param               = topoCoeff{m};
%         sweT                = topo_full.(glacier).(param)*mlrCoeff(m);
%         sweBMS.(glacier)    = sweBMS.(glacier) + sweT;
%     end
%         %Set min to 0
%     sweBMS.(glacier)(sweBMS.(glacier)<0) = 0;
% 
%     %RMSE
% %     sampledtemp = sweBMS.(glacier)(options.ENgrid.(glacier)(:,2),options.ENgrid.(glacier)(:,1));
% %     estGrid     = diag(sampledtemp);
% 
%         %DataObs_RMSEeven.(type).(glacier)(n,mc) = sqrt(mean((estGrid-realGrid.(glacier)(:,1)).^2));
% end
% 
% sweBMS_alldata = sweBMS;
% save('PaperII_AblationArea.mat','sweBMS_alldata', '-append')
%% Blank out accum area %

load PaperII_AblationArea.mat sweBMS_alldata

% Glacier 4
AblationArea.G4 = sweBMS_alldata.G4;
x = repmat(1:options.mapsize(1,2),options.mapsize(1,1),1);
y = repmat((1:options.mapsize(1,1))',1,options.mapsize(1,2));

ELA.G4 = [60,19;54,20;50,21;35,23;28,29;29,36;30,38];
for i = 2:length(ELA.G4)
    m = (ELA.G4(i,1)-ELA.G4(i-1,1))/(ELA.G4(i,2)-ELA.G4(i-1,2));
    b = ELA.G4(i,1)-m*ELA.G4(i,2);

    y_test = m*x+b;
    y_nan = y_test>y;
    AblationArea.G4(y_nan) = NaN;
end

T = AblationArea.G4 ~= sweBMS_alldata.G4;
AblationArea.G4(T) = -0.1;
AblationArea.G4(options.mapNaN.G4) = NaN;

% clf; imagesc(AblationArea.G4); hold on
% plot(ELA.G4(:,2),ELA.G4(:,1),'k--')


% Glacier 2
AblationArea.G2 = sweBMS_alldata.G2;
x = repmat(1:options.mapsize(2,2),options.mapsize(2,1),1);
y = repmat((1:options.mapsize(2,1))',1,options.mapsize(2,2));

ELA.G2 = [90,67;87,80;82,89;74,90;68,90];
for i = 2:length(ELA.G2)
    m = (ELA.G2(i,1)-ELA.G2(i-1,1))/(ELA.G2(i,2)-ELA.G2(i-1,2));
    b = ELA.G2(i,1)-m*ELA.G2(i,2);

    y_test = m*x+b;
    y_nan = y_test<y;
    AblationArea.G2(y_nan) = NaN;
end
T = AblationArea.G2 ~= sweBMS_alldata.G2;
AblationArea.G2(T) = -0.1;
AblationArea.G2(options.mapNaN.G2) = NaN;

% clf; imagesc(AblationArea.G2); hold on
% plot(ELA.G2(:,2),ELA.G2(:,1),'k--')

% Glacier 13
AblationArea.G13 = sweBMS_alldata.G13;
x = repmat(1:options.mapsize(3,2),options.mapsize(3,1),1);
y = repmat((1:options.mapsize(3,1))',1,options.mapsize(3,2));

ELA.G13 = [95,61;105,74;99,89;87,90];
for i = 2:length(ELA.G13)
    m = (ELA.G13(i,1)-ELA.G13(i-1,1))/(ELA.G13(i,2)-ELA.G13(i-1,2));
    b = ELA.G13(i,1)-m*ELA.G13(i,2);

    y_test = m*x+b;
    y_nan = y_test<y;
    AblationArea.G13(y_nan) = NaN;
end

extra_ELA = [60,28;73,46];
i = 2;
m = (extra_ELA(i,1)-extra_ELA(i-1,1))/(extra_ELA(i,2)-extra_ELA(i-1,2));
b = extra_ELA(i,1)-m*extra_ELA(i,2);
y_test = m*x+b;
y_nan(:,24:45) = y_test(:,24:45)<y(:,24:45);
AblationArea.G13(y_nan) = NaN;

T = AblationArea.G13 ~= sweBMS_alldata.G13;
AblationArea.G13(T) = -0.1;
AblationArea.G13(options.mapNaN.G13) = NaN;

% clf; imagesc(AblationArea.G13); hold on
% plot(ELA.G13(:,2),ELA.G13(:,1),'k--')

PlotTopoParameter(AblationArea, 'name', 'B_W (m w.e.)',fullSWE.S2.input, 'black','massB')





