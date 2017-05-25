%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compare Signals (POI1) of Lungsim B-files with Spiroware B-files
%
% Version 4.0, 17. Sept. 2013
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=compareSignals(signalLungsim,signalSpiroware,figureNo,parameters)

    spirowareType = parameters.Operation.spirowareType;
    if spirowareType == 0
        nGraphs=4;
    else
        nGraphs=2;
    end    

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Shift zero of time scale
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    startLungsim=find(signalLungsim.O2Poi1>0.9,1);
    startSpiroware=find(signalSpiroware.O2>0.9,1);
    timeOffset=signalLungsim.ts(startLungsim)-signalSpiroware.ts(startSpiroware);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graphical output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(10000+figureNo)
    hSubplot10000(1)=subplot(nGraphs,1,1);
    plot(signalLungsim.ts,signalLungsim.O2Poi1,'b',signalSpiroware.ts+timeOffset,signalSpiroware.O2,'r','Parent',hSubplot10000(1));
    title('O2 signals!')
    xlabel('t / s')
    ylabel('O2 / -')
    legend('Lungsim','Spiroware')

    hSubplot10000(2)=subplot(nGraphs,1,2);
    plot(signalLungsim.ts,signalLungsim.CO2Poi1,'b',signalSpiroware.ts+timeOffset,signalSpiroware.CO2,'r','Parent',hSubplot10000(2));
    title('CO2 signals!')
    xlabel('t / s')
    ylabel('CO2 / -')
    legend('Lungsim','Spiroware')

    if spirowareType == 0
        hSubplot10000(3)=subplot(nGraphs,1,3);
        plot(signalLungsim.ts,signalLungsim.N2Poi1,'b',signalSpiroware.ts+timeOffset,signalSpiroware.N2,'r','Parent',hSubplot10000(3));
        title('N2 signals!')
        xlabel('t / s')
        ylabel('N2 / -')
        legend('Lungsim','Spiroware')
        
        hSubplot10000(4)=subplot(nGraphs,1,4);
        plot(signalLungsim.ts,signalLungsim.N2Poi1+signalLungsim.O2Poi1+signalLungsim.CO2Poi1,'b',signalSpiroware.ts+timeOffset,signalSpiroware.N2+signalSpiroware.O2+signalSpiroware.CO2,'r','Parent',hSubplot10000(4));
        title('Sum signals!')
        xlabel('t / s')
        ylabel('N2+O2+CO2 / -')
        legend('Lungsim','Spiroware')
        
    end
    linkaxes(hSubplot10000,'x')
    
    figure(20000+figureNo)
    plot(signalLungsim.ts,signalLungsim.IvFilter,'b',signalSpiroware.ts+timeOffset,signalSpiroware.Iv,'r');
    title('Flow signals!')
    xlabel('t / s')
    ylabel('Flow / m^3/s')
    legend('Lungsim','Spiroware')
end