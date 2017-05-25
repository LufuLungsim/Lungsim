%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Selecting interval of input file
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.0, 27. März 2008
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tmin,tmax]=cropDataSpirowareMassSpec(signal,parameters,header)
    fprintf('\nmanual calibration\n')
    mmAir           =   parameters.Operation.MMair*0.95;
    mm2Percent      =   parameters.Operation.MM2Percent;
    
    torontoFile     =   parameters.Simulation.torontoFile;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graphical output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hfig=figure(1000);
    set(hfig,'Name',header);
    set(hfig,'NumberTitle','off');
    hold off
    if torontoFile
        plot(signal.ts,signal.CO2*100,'b',signal.ts,signal.He*100,'g',signal.ts,signal.SF6*100,'r-');
        legendStrings={'CO2','O2','MMss (arb units)'};
        axis([-inf,inf,-1,9])
    else
        plot(signal.ts,signal.CO2*100,'b',signal.ts,signal.O2*100,'g',signal.ts,(signal.MMss-mmAir)/mm2Percent*100,'r-');
        legendStrings={'CO2','O2','MMss (arb units)'};
        axis([-inf,inf,-10,100])
    end
    title('left/right mouse button to chose interval, <enter> to proceed)')
    legend(legendStrings)
    xlabel('t / s')
    ylabel('x / %')
    
    fprintf(['indicate lower/upper bound of ',header,' with left/right mouse button; press <enter> to confirm.\n'])
    
    tmin=signal.ts(1);
    tmax=signal.ts(end);
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
                    hmin=plot([tmin,tmin],[-10,40],'m--');
                case 3
                    tmax=t;
                    if hmax~=0
                        delete(hmax);
                    end
                    hmax=plot([tmax,tmax],[-10,40],'m-.');
                otherwise
                    fprintf('unknown button/key pressed\n')
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
    
    
