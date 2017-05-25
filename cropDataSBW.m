%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Selecting interval of input file
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.1, 27. August 2012
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmin,tmax]=cropDataSBW(signal,breathIndex,breathTimes,parameters,header)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Data setup
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minIndex = breathIndex(1);
    maxIndex = breathIndex(end);
    ts = signal.ts(minIndex:maxIndex);
    MMss = signal.MMssPoi3(minIndex:maxIndex);
    MMssCalc = signal.MMssCalc(minIndex:maxIndex);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graphical output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hfig=figure(1000);
    set(hfig,'Name',header);
    set(hfig,'NumberTitle','off');
    hold off
    plot(ts,MMss*1000,'b',ts,MMssCalc*1000,'g');
    hold on
    minMM=min([MMss;MMssCalc])*1000;
    maxMM=max([MMss;MMssCalc])*1000;
    for time=breathTimes
        plot([time,time],[minMM,maxMM],'k-.');
    end
    hold off
    title('left/right mouse button to chose SBW analysis interval, <enter> to proceed)')
    legendStrings={'MMss','MMssCalc'};
    legend(legendStrings)
    xlabel('t / s')
    ylabel('x / %')
    axis([-inf,inf,minMM,maxMM])
    
    
    disp(['indicate lower/upper bound of ',header,' with left/right mouse button; press <enter> to confirm.'])
    
    tmin=ts(1);
    tmax=ts(end);
    hmin=0;
    hmax=0;
    button=1;
    while(button)
        [t,y,button]=ginput(1);
        hold on
        if ~isempty(button)
            switch button
                case 1
                    tmin=t;
                    if hmin~=0
                        delete(hmin);
                    end
                    hmin=plot([tmin,tmin],[minMM,maxMM],'m--');
                case 3
                    tmax=t;
                    if hmax~=0
                        delete(hmax);
                    end
                    hmax=plot([tmax,tmax],[minMM,maxMM],'m-.');
                otherwise
                    disp('unknown button/key pressed')
            end
        end
        hold off
    end
    
    if tmin>tmax
        tdummy=tmax;
        tmax=tmin;
        tmin=tdummy;
    end
    
    close(hfig);
end
    
    
