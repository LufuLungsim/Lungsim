%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to correct signals in various modes
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.4, 20. Oktober 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal,breathIndex]=signalCorrectionWbreath(signal,parameters,flag)

    dataType	=   parameters.Simulation.SF6;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Breath segmentation and calculation of volume (breath detected variant)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [breathIndex,breathTimes]=breathDetectionSimple(signal.ts,signal.Vol);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag~=0
        
        hfig=figure(10+flag);
        set(hfig,'Name','Signals');
        hSubplot10(1)=subplot(2,1,1);
        plot(signal.ts,signal.Iv*1000,'b',signal.ts,signal.Vol*1000,'g','Parent',hSubplot10(1));
        hold on
        minVol = (floor((min(min(signal.Iv),min(signal.Vol))*1000))*10)/10;
        maxVol = -(floor(-(max(max(signal.Iv),max(signal.Vol))*1000))*10)/10;
        for time=breathTimes
            plot([time,time],[minVol,maxVol],'k-.');
        end
        hold off
        title('raw data')
        xlabel('t / s')
        ylabel('Iv,V / l/s,l')
        legend('Iv', 'Vol')
        axis([-inf,inf,-inf,inf])
        
        hSubplot10(2)=subplot(2,1,2);
        plot(signal.ts,signal.MM*1000,'b','Parent',hSubplot10(2));
        if dataType == 0
            hold on
            plot(signal.ts,signal.MMss*1000,'r','Parent',hSubplot10(2));
            hold off
        end
        title('MM signal')
        xlabel('t / s')
        ylabel('MM / g/mol')
        if dataType == 0
            legend('MM','MMss')
            MMmin = 22;
            MMmax = 30;
        else
            legend('MM')
            MMmin = 28;
            MMmax = 36;
        end
        hold on
        for time=breathTimes
            plot([time,time],[MMmin,MMmax],'k-.');
        end
        hold off
        axis([-inf,inf,-inf,inf])
        linkaxes(hSubplot10,'x')
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function to detect breaths and its indices, works only for wbreath that
% contains breath-wise information!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [index,breaths]=breathDetectionSimple(time,vol)

    signChange=((sign(vol(1:end-1))<1).*(sign(vol(2:end))==1))+((sign(vol(1:end-1))>-1).*(sign(vol(2:end))==-1));
    index=find(signChange)+1;

    breaths = time(index);                               % determine breath times
end



