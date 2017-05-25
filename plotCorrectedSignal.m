%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to plot toronto signals 
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 4.8, 1. December 2015
% Original by Markus Roos, NM GmbH, adapted by Jerry Wolfensberger
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  plotCorrectedSignal( signal, parameters, breathTimes, flag )

    torontoFile = parameters.Simulation.torontoFile;
    
    if torontoFile
        hfig=figure(10+flag);
        set(hfig,'Name','Signals');
        hSubplot10(1)=subplot(2,1,1);
        plot(signal.ts,signal.IvFilter, signal.ts,signal.VolFilter,...
            signal.ts,signal.IvEffBTPS,signal.ts,signal.VolEffBTPS,signal.ts,signal.VolBTEffBTPS,...
            signal.ts,signal.delayLineFilter*1e-3,'k--','Parent',hSubplot10(1));
        hold on
        for time=breathTimes
            plot([time,time],[-1.5e-3,1.5e-3],'k-.');
        end
        hold off
        title('delayed data for air flow')
        xlabel('t / s')
        ylabel('Iv, V / m^3/s, m^3')
        legend('IvFilter', 'VolFilter', 'IvEffBTPS', 'VolEffBTPS', 'VolBTEffBTPS', 'delayLine')
        axis([-inf,inf,-inf,inf])
        
        hSubplot10(2)=subplot(2,1,2);
        plot(signal.ts,signal.He,'b',signal.ts,signal.SF6,'r',signal.ts,signal.CO2,'g','Parent',hSubplot10(2));
        title('He, SF6, and CO2 signals')
        xlabel('t / s')
        ylabel('x / -')
        legend('He', 'SF6', 'CO2')
        hold on
        for time=breathTimes
            plot([time,time],[0,0.1],'k-.');
        end
        hold off
        axis([-inf,inf,-inf,inf])
        linkaxes(hSubplot10,'x')
        
    else    % normal Spiroware based file
        hfig=figure(10+flag);
        set(hfig,'Name','Corrected signals');
        
        hSubplot10(1)=subplot(3,1,1);
        plot(signal.ts,signal.IvDelay, signal.ts,signal.VolFilter,...
            signal.ts,signal.IvEffBTPS,signal.ts,signal.VolEffBTPS,signal.ts,signal.VolBTEffBTPS,...
            'k--','Parent',hSubplot10(1));
        hold on
        for time=breathTimes
            plot([time,time],[-1.5e-3,1.5e-3],'k-.');
        end
        hold off
        title('delayed data for air flow')
        xlabel('t / s')
        ylabel('Iv, V / m^3/s, m^3')
        legend('IvFilter', 'VolFilter', 'IvEffBTPS', 'VolEffBTPS', 'VolBTEffBTPS', 'delayLine')
        axis([-inf,inf,-inf,inf])
        
        hSubplot10(2)=subplot(3,1,2);
        plot(signal.ts,signal.MMssPoi3*1000,'Parent',hSubplot10(2));
        hold on
        for time=breathTimes
            plot([time,time],[28,33],'k-.');
        end
        hold off
        title('MMss signal at POI3')
        xlabel('t / s')
        ylabel('MM / g/mol')
        legend('MMssPOI3')
        axis([-inf,inf,-inf,inf])
        
        hSubplot10(3)=subplot(3,1,3);
        plot(signal.ts,signal.CO2Poi1, signal.ts,signal.O2Poi1,'Parent',hSubplot10(3));
        title('CO2 and O2 signals at POI3')
        xlabel('t / s')
        ylabel('x / -')
        legend('CO2POI1', 'O2POI1')
        hold on
        for time=breathTimes
            plot([time,time],[0,1.0],'k-.');
        end
        hold off
        axis([-inf,inf,-inf,inf])
        %        pause(30)
        linkaxes(hSubplot10,'x')
        
        if parameters.Simulation.verb
            figure(20+flag);
            subplot(2,1,1);
            plot(signal.ts,signal.CO2Poi3, signal.ts,signal.O2Poi3, signal.ts,signal.N2Poi3, signal.ts,signal.ArPoi3);
            hold on
            for time=breathTimes
                plot([time,time],[0.0,1.0],'k-.');
            end
            hold off
            title('molar fractions at POI3')
            xlabel('t / s')
            ylabel('x / -')
            legend('CO2', 'O2', 'N2', 'Ar')
            axis([-inf,inf,-inf,inf])
            
            subplot(2,1,2);
            plot(signal.ts,signal.MMssPoi3*1000, signal.ts,signal.MMssCalc*1000);
            hold on
            for time=breathTimes
                plot([time,time],[28,33],'k-.');
            end
            hold off
            title('measured (filtered) and calculated MMss at POI3')
            xlabel('t / s')
            ylabel('MM / g/mol')
            legend('MMssPOI3', 'MMssCalc')
            axis([-inf,inf,-inf,inf])
            
            figure(40+flag);
            subplot(1,1,1);
            %    plot(signal.ts,signal.IvssFilter*1000,signal.ts,signal.IvssSmoothFilter*1000,signal.ts,signal.IvssSmooth1Filter*1000);
            plot(signal.ts,signal.IvssSmoothFilter*1000,signal.ts,signal.IvssSmooth1Filter*1000);
            hold on
            for time=breathTimes
                plot([time,time],[-4.0e-3,0.0],'k-.');
            end
            hold off
            title('Variation of Ivss')
            xlabel('t / s')
            ylabel('Ivss / -')
            %    legend('raw', 'smooth', 'moving average')
            legend('smooth', 'moving average')
            axis([-inf,inf,-inf,inf])
        end
        
        hfig=figure(30+flag);
        
        figureSize=get(hfig,'Position');    % rescaling graph if necessary (for better visibility)
        height=figureSize(4);
        width=figureSize(3);
        if width>height
            figureSize(4)=1.6*height;
            figureSize(2)=figureSize(2)-0.6*height;
            set(hfig,'Position',figureSize);
        end
        
        set(hfig,'Name','Species concentration and flow signals at POI1');
        hSubplot30(1)=subplot(2,1,1);
        plot(signal.ts,signal.CO2Poi1*100, signal.ts,signal.O2Poi1*100, signal.ts,signal.N2Poi1*100,...
            signal.ts,signal.ArPoi1*100,'Parent',hSubplot30(1));
        hold on
        for time=breathTimes
            plot([time,time],[0.0,100.0],'k-.');
        end
        hold off
        title('molar fractions at POI1')
        xlabel('t / s')
        ylabel('x / %')
        legend('CO2', 'O2', 'N2', 'Ar')
        axis([-inf,inf,-inf,inf])
        
        hSubplot30(2)=subplot(2,1,2);
        plot(signal.ts,signal.IvEffBTPS*1e6,signal.ts,signal.VolBTEffBTPS*1e6,'Parent',hSubplot30(2));
        hold on
        for time=breathTimes
            plot([time,time],[-1.5e3,1.5e3],'k-.');
        end
        hold off
        title('flow related data')
        xlabel('t / s')
        ylabel('I_V_,_B_T_P_S, V_B_T_P_S / ml/s, ml')
        legend('IvEffBTPS', 'VolBTEffBTPS')
        axis([-inf,inf,-inf,inf])
        linkaxes(hSubplot30,'x')
        
        
        hfig=figure(50+flag);
        set(hfig,'Name','Comparison of  MMss and MMssCalc (POI3)');
        hSubplot50(1)=subplot(2,1,1);
        plot(signal.ts,signal.MMssPoi3*1000,'b',signal.ts,signal.MMssCalc*1000,'r','Parent',hSubplot50(1));
        hold on
        minMM=min([signal.MMssPoi3;signal.MMssCalc])*1000;
        maxMM=max([signal.MMssPoi3;signal.MMssCalc])*1000;
        for time=breathTimes
            plot([time,time],[minMM,maxMM],'k-.');
        end
        hold off
        title('Molar masses')
        xlabel('t / s')
        ylabel('MMss / g/mol')
        legendStrings={'MMss', 'MMssCalc'};
        legend(legendStrings)
        axis([-inf,inf,-inf,inf])
        
        hSubplot50(2)=subplot(2,1,2);
        plot(signal.ts,signal.DiffMMss*1000,'b',signal.ts,zeros(size(signal.ts)),'k','Parent',hSubplot50(2));
        minMMDiff=min(signal.DiffMMss)*1000;
        maxMMDiff=max(signal.DiffMMss)*1000;
        hold on
        for time=breathTimes
            plot([time,time],[minMMDiff,maxMMDiff],'k-.');
        end
        hold off
        title('Molar mass difference')
        xlabel('t / s')
        ylabel('DiffMMss / g/mol')
        legend('DiffMMss')
        axis([-inf,inf,-inf,inf])
        linkaxes(hSubplot50,'x')
        
        figure(30+flag); % show most important info (at this point of operation
    end
    
end

