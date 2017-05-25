%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to transform signal data to Leicester format
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.13, 10. Mai 2012
% Markus Roos, LuFu, Inselspital Bern
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [leicesterSignal] = outputLeicester(signal,verb)

    nData = floor(length(signal.ts)/2);
    
    flow	=   signal.IvEffBTPS;
    vol 	=   signal.VolEffBTPS;
    volBT	=   signal.VolBTEffBTPS;
    
    minLeicester=zeros(nData,1);
    secLeicester=zeros(nData,1);
    milliLeicester=zeros(nData,1);
    flowLeicester=zeros(nData,1);
    volumeLeicester=zeros(nData,1);
    volumeBTLeicester=zeros(nData,1);
    n2Leicester=zeros(nData,1);
    ts=zeros(nData,1);
    
    tsZero=(signal.ts(1)+signal.ts(2))/2;
    for i = 1:nData
        ts(i)=(signal.ts(2*i-1)+signal.ts(2*i))/2-tsZero;
        [minutes,seconds,millisec]=getTimeFormat(ts(i));
        minLeicester(i)=minutes;
        secLeicester(i)=seconds;
        milliLeicester(i)=millisec;
        flowLeicester(i)=-1000*(flow(2*i-1)+flow(2*i))/2;
        volumeLeicester(i)=-1000*(vol(2*i-1)+vol(2*i))/2;
        volumeBTLeicester(i)=-1000*(volBT(2*i-1)+volBT(2*i))/2;
        n2Leicester(i)=100*(signal.N2Poi1(2*i-1)+signal.N2Poi1(2*i))/2;
    end
    
    leicesterSignal.min=minLeicester;
    leicesterSignal.sec=secLeicester;
    leicesterSignal.milli=milliLeicester;
    leicesterSignal.flow=flowLeicester;
    leicesterSignal.volume=volumeLeicester;
    leicesterSignal.volumeBT=volumeBTLeicester;
    leicesterSignal.N2=n2Leicester;
    
    if verb
        i=0;
        startIndex = 1;
        actualSign = sign(volumeBTLeicester(1));
        colorCode = 'bkr';
        while i<length(volumeBTLeicester)
            i=i+1;
            newSign = sign(volumeBTLeicester(i));
            if newSign ~= actualSign
                endIndex = i-1;
                figure(900);
                plot(ts(startIndex:endIndex),volumeLeicester(startIndex:endIndex),colorCode(actualSign+2));
                hold on;
                figure(901);
                if actualSign == -1
                    plot(volumeBTLeicester(startIndex:endIndex),n2Leicester(startIndex:endIndex),'b.');
                else
                    maxVolume = max(volumeBTLeicester(startIndex:endIndex));
                    plot(volumeBTLeicester(startIndex:endIndex)-maxVolume,n2Leicester(startIndex:endIndex),'r.');
                end
                hold on;
                actualSign = newSign;
                startIndex = endIndex+1;
            end
        end
        figure(900);
        plot(ts,n2Leicester/10,'k');
        title('Leicester-like graph')
        xlabel('t / s')
        ylabel('vol / l, N2 / %/10')
        axis([-inf,inf,-1,8]);
        hold off;
        
        figure(901);
        title('Breath compilation')
        xlabel('vol / l')
        ylabel('N2 / %')
        axis([-inf,inf,0,80]);
        hold off;
    end
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time format conversion
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [minutes,seconds,millisec]=getTimeFormat(t)

time=floor(t);
millisec=round(1000*(t-time));
seconds=mod(time,60);
minutes=(time-seconds)/60;

end

