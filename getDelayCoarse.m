%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to filter signals
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.0, 25. August 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [delayO2,delayMMss] = getDelayCoarse(signal,parameters)

    dt      =	parameters.Simulation.dt;
    verb    =   parameters.Simulation.verb;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % triming of signals (to make sure that only signal portions are considered
    % that enable for all signals (i.e., the also the latest MMss, etc.)
    % ------ to be implemented -------------------
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % finding instant of steepest ascent/descent in O2, CO2 and MMss
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    O2=signal.O2;
    O2=-(O2-mean(O2));
    O2=O2/max(abs(O2));
    MMss=signal.MMss;
    MMss=MMss-mean(MMss);
    MMss=MMss/max(abs(MMss));
    CO2=signal.CO2;
    CO2=CO2-mean(CO2);
    CO2=CO2/max(abs(CO2));
    
    xCorrelationCO2MMss         =   real(ifft(fft(MMss).*conj(fft(CO2))));
    [maxCO2MMss,indexCO2MMss]   =   max(xCorrelationCO2MMss);
    xCorrelationCO2O2           =   real(ifft(fft(O2).*conj(fft(CO2))));
    [maxCO2O2,indexCO2O2]       =   max(xCorrelationCO2O2);
    
    if indexCO2O2<length(signal.ts)/2
        delayO2 = indexCO2O2*0.005;
    else
        delayO2 = (indexCO2O2-length(signal.ts))*dt;
    end
    
    if indexCO2MMss<length(signal.ts)/2
        delayMMss = indexCO2MMss*dt;
    else
        delayMMss = (indexCO2MMss-length(signal.ts))*dt;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optional output to control algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if verb>0
        figure(401)
        subplot(2,1,1)
        plot(signal.ts,[O2,CO2,MMss]);
        legend('O2scaled','CO2scaled','MMssscaled');
        subplot(2,1,2)
        plot([xCorrelationCO2O2,xCorrelationCO2MMss]);
        legend('Corr CO2-O2','Corr CO2-MMss');
    end
end
