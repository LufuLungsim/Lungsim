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
function [tmin,tmax]=cropDataWbreath(signal,parameters,header)

    dataType        =   parameters.Simulation.SF6;
    
    if dataType == 0
        MMmin = 22;
        MMmax = 30;
    else
        MMmin = 28;
        MMmax = 36;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graphical output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hfig=figure(1000);
    set(hfig,'Name',header);
    set(hfig,'NumberTitle','off');
    hold off
    plot(signal.ts,signal.MM*1000,'b');
    legendStrings={'MM'};
    axis([-inf,inf,MMmin,MMmax])
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
                    hmin=plot([tmin,tmin],[MMmin,MMmax],'m--');
                case 3
                    tmax=t;
                    if hmax~=0
                        delete(hmax);
                    end
                    hmax=plot([tmax,tmax],[MMmin,MMmax],'m-.');
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


