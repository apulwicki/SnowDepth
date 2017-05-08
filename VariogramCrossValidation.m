%% Variogram - transect
 %%Fitting function cannot have small values for the variance!
runs = 10;%00;
multFactor = 10000;
figure(1); clf

for g = 1:3    
    glacier = options.glacier{g};
    range.(glacier) = zeros(runs,3); weights.(glacier) = range.(glacier); 
    sill.(glacier) = range.(glacier); nugget.(glacier) = range.(glacier);
    
    
for i = 1:runs
Iin     = randi([1 length(SWE(g).swe)],round(length(SWE(g).swe)*2/3),1);
Iout    = setdiff( 1:length(SWE(g).swe), Iin)';
 
%Calculate variogram parameters for IN data
    variogramIN.(glacier) = variogramAlex([SWE(g).swe(Iin) SWE(g).utm(Iin,1:2)], 15, 'default');
        Ivario = variogramIN.(glacier).num~=0;
        variogramIN.(glacier).val         = variogramIN.(glacier).val(Ivario)*multFactor;
        variogramIN.(glacier).meanDist    = variogramIN.(glacier).meanDist(Ivario);
        variogramIN.(glacier).num         = variogramIN.(glacier).num(Ivario);
        variogramIN.(glacier).binCentre   = variogramIN.(glacier).binCentre(Ivario);
    fit.(glacier) = variofitAlex(variogramIN.(glacier),glacier,6);
        range.(glacier)(i,:)      = [fit.(glacier).spherical.range, fit.(glacier).exponential.range, fit.(glacier).gaussian.range];
        sill.(glacier)(i,:)       = [fit.(glacier).spherical.sill, fit.(glacier).exponential.sill, fit.(glacier).gaussian.sill];
        nugget.(glacier)(i,:)     = [fit.(glacier).spherical.nugget, fit.(glacier).exponential.nugget, fit.(glacier).gaussian.nugget];

                Imal = range.(glacier)>5000 & range.(glacier)<80;
            range.(glacier)(Imal) = NaN;    nugget.(glacier)(Imal) = NaN; 
            sill.(glacier)(Imal) = NaN;     weights.(glacier)(Imal) = NaN;     

%Calculate variogram parameters for OUT data
    variogramOUT.(glacier) = variogramAlex([SWE(g).swe(Iout) SWE(g).utm(Iout,1:2)], 15, 'default');
    variogramOUT.(glacier).val = variogramOUT.(glacier).val*multFactor;
    
    %Prediction
    lags.(glacier) = variogramOUT.(glacier).meanDist;
     %Spherical
        pred.(glacier).spherical        = nugget.(glacier)(i,1) + ( sill.(glacier)(i,1)*( 1.5*(lags.(glacier)/range.(glacier)(i,1)) - 0.5*(lags.(glacier)/range.(glacier)(i,1)).^3)); 
        pred.(glacier).spherical        = [pred.(glacier).spherical(lags.(glacier)<range.(glacier)(i,1)); (sill.(glacier)(i,1)+nugget.(glacier)(i,1))*(lags.(glacier)>range.(glacier)(i,1))];
        pred.(glacier).spherical(pred.(glacier).spherical == 0) = [];
     %Exponential
        pred.(glacier).exponential      = nugget.(glacier)(i,2) + sill.(glacier)(i,2)*(1-exp(-lags.(glacier)/range.(glacier)(i,2))); 
     %Exponential
        pred.(glacier).gaussian         = nugget.(glacier)(i,3) + ( sill.(glacier)(i,3)*( 1 - exp(-3*lags.(glacier).^2/range.(glacier)(i,3)^2))); 
    
    %RMSE
    for t = 1:3
        if      t==1; type = 'spherical';
        elseif  t==2; type = 'exponential';
        elseif  t==3; type = 'gaussian';
        end
        varioRMSE.(glacier)(i,t) = sqrt(mean((pred.(glacier).(type)(:,1)-variogramOUT.(glacier).val).^2));
    end
end

%Weights for params
    weightRMSE.(glacier) = 1./varioRMSE.(glacier);
    weightRMSE.(glacier) = weightRMSE.(glacier)./repmat(sum(weightRMSE.(glacier)),runs,1);
    
%Mean with weights
    Mrange.(glacier)  = nansum(range.(glacier).*weightRMSE.(glacier));
    Mnugget.(glacier) = nansum(nugget.(glacier).*weightRMSE.(glacier));
    Msill.(glacier)   = nansum(sill.(glacier).*weightRMSE.(glacier));

%Prediction with final parameters
    variogramALL.(glacier) = variogramAlex([SWE(g).swe SWE(g).utm(:,1:2)], 15, 'default');
    lags.(glacier) = variogramALL.(glacier).meanDist;
     %Spherical
        pred.(glacier).spherical        = Mnugget.(glacier)(1) + ( Msill.(glacier)(1)*( 1.5*(lags.(glacier)/Mrange.(glacier)(1)) - 0.5*(lags.(glacier)/Mrange.(glacier)(1)).^3)); 
        pred.(glacier).spherical        = [pred.(glacier).spherical(lags.(glacier)<Mrange.(glacier)(1)); (Msill.(glacier)(1)+Mnugget.(glacier)(1))*(lags.(glacier)>Mrange.(glacier)(1))];
        pred.(glacier).spherical(pred.(glacier).spherical == 0) = [];
     %Exponential
        pred.(glacier).exponential      = Mnugget.(glacier)(2) + Msill.(glacier)(2)*(1-exp(-lags.(glacier)/Mrange.(glacier)(2))); 
     %Exponential
        pred.(glacier).gaussian         = Mnugget.(glacier)(3) + ( Msill.(glacier)(3)*( 1 - exp(-3*lags.(glacier).^2/Mrange.(glacier)(3)^2))); 
    
subplot(1,3,g)
plot(lags.(glacier), pred.(glacier).spherical); hold on
plot(lags.(glacier), pred.(glacier).exponential);
plot(lags.(glacier), pred.(glacier).gaussian);
plot(variogramOUT.(glacier).meanDist, variogramOUT.(glacier).val,'o')

%Final weights for the three fits
    for t = 1:3
        if      t==1; type = 'spherical';
        elseif  t==2; type = 'exponential';
        elseif  t==3; type = 'gaussian';
        end
        weightRMSE.final(g,t) = sqrt(mean((pred.(glacier).(type)(:,1)-variogramALL.(glacier).val).^2));
    end
    weightRMSE.final = 1./weightRMSE.final;
    weightRMSE.final = weightRMSE.final./repmat(sum(weightRMSE.final,2),1,3);
    
    pred.(glacier).final = pred.(glacier).spherical*weightRMSE.final(g,1)+...
                           pred.(glacier).exponential*weightRMSE.final(g,2)+...
                           pred.(glacier).gaussian*weightRMSE.final(g,3);
    
plot(lags.(glacier), pred.(glacier).final);
    legend('Spherical','Exponential','Gaussian','Validation Data','Final','Location','best')
  
end
    %saveFIG(['/home/glaciology1/Documents/Data/Plots/variofull',glacier])
