%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function for plot only modus to check for maligne data sets
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 14. Mai 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=plotSignalWbreath(signal,parameters,flag)

dataType	=   parameters.Simulation.SF6;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Plotting results
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if flag~=0
    hfig=figure(10+flag);
    set(hfig,'Name','Raw signals');
    
    hSubplot10(1)=subplot(2,1,1);
    plot(signal.ts,signal.Iv*1000,signal.ts,signal.Vol*1000);
    xlabel('t / s')
    ylabel('Iv / l/s, Vol / l')
    legend('Iv','Vol')
    axis([-inf,inf,-inf,inf])
    
    hSubplot10(2)=subplot(2,1,2);
    plot(signal.ts,signal.MM*1000,'b');
    if dataType == 0 
        hold on
        plot(signal.ts,signal.MMss*1000,'r');
        hold off
    end
    title('Molar mass signal')
    xlabel('t / s')
    ylabel('MM / g/mol')
    if dataType == 0 
        legend('MM','MMss')
    else
        legend('MM')
    end
    axis([-inf,inf,-inf,inf])
    
    linkaxes(hSubplot10,'x')
    
end
