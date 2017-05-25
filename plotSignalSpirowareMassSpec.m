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
function []=plotSignalSpirowareMassSpec(signal,parameters,flag)

    torontoFile     =   parameters.Simulation.torontoFile;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag~=0
        if ~torontoFile
            hfig=figure(10+flag);
            set(hfig,'Name','Raw signals');
            
            subplot(3,1,1);
            plot(signal.ts,signal.Iv*1000);
            xlabel('t / s')
            ylabel('Iv / l/s')
            legend('Iv')
            axis([-inf,inf,-inf,inf])
            
            subplot(3,1,2);
            plot(signal.ts,signal.CO2*100,'g',signal.ts,signal.O2*100,'b');
            title('Gas species')
            xlabel('t / s')
            ylabel('x / %')
            legend('CO2','O2')
            axis([-inf,inf,-inf,inf])
            
            subplot(3,1,3);
            plot(signal.ts,signal.MMss*1000,'g',signal.ts,signal.MMms*1000,'b');
            title('Molar mass signal')
            xlabel('t / s')
            ylabel('x / g/mol')
            legend('MMss','MMms')
            axis([-inf,inf,-inf,inf])
            
        else
            hfig=figure(10+flag);
            set(hfig,'Name','Raw Signals');
            subplot(2,1,1);
            plot(signal.ts,signal.Iv*1000);
            xlabel('t / s')
            ylabel('Iv / l/s')
            legend('Iv')
            axis([-inf,inf,-inf,inf])
            
            subplot(2,1,2);
            plot(signal.ts,signal.He*100,'b',signal.ts,signal.SF6*100,'r',signal.ts,signal.CO2*100,'g');
            title('He, SF6, and CO2 signals')
            xlabel('t / s')
            ylabel('x / %')
            legend('He', 'SF6', 'CO2')
            axis([-inf,inf,-inf,inf])
        end
    end
end
