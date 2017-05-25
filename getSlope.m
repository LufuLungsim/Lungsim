%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interactive calculation of slope, works also for batch operation
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 14. Mai 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fitParameters,interactive,error]=getSlope(type,quantityString,minRelative,maxRelative,x,y,flow,interactive,rmsCrit)

    valueMax=max(x);
    valueMin=min(x);
    
    if minRelative==0
        xMin=valueMin;
    else
        xMin=minRelative*valueMax;
    end
    if maxRelative==0
        xMax=valueMax;
    else
        xMax=maxRelative*valueMax;
    end
    
    if strcmp(type,'normlin')   % normalization of y-data if necessary
        indices=find(1==(x>xMin).*(x<xMax));
        y=y.*(1/mean(y(indices)));
    end
    
    hfig=figure(901);
    if interactive>1    % adjusting figure size if interactive
        figureSize=get(hfig,'Position');
        height=figureSize(4);
        figureSize(4)=1.6*height;
        figureSize(2)=figureSize(2)-0.6*height;
        set(hfig,'Position',figureSize);
    end
    set(hfig,'Name',strcat('Slope determination (method: ',type,')'));
    set(hfig,'NumberTitle','off');
    
    subplot(2,1,2)
    plot(x,flow,'b');
    title('volume flow plot');
    xlabel('Vol / m^3');
    ylabel('Iv / m^3/');
    axis([valueMin,valueMax,-inf,inf]);
    
    subplot(2,1,1)
    hold off
    plot(x,y,'b');
    title('Indicate lower/upper bound; <enter> confirms, <ESC> batch mode.');
    xlabel('Vol / m^3');
    ylabel(quantityString);
    hold on
    
    yMin=min(y);
    yMax=max(y);
    yDiff=yMax-yMin;
    yMin=yMin-yDiff/2;
    yMax=yMax+yDiff/2;

    [coeff,fit,rms,rSquared,error]=fitSlope(type,xMin,xMax,x,y);
    if error
%        fitParameters=coeff;
        fitParameters=zeros(1,2);
        close(hfig);
        return;
    end

    hfit=plot(x,fit,'r-.');
    hmin=plot([xMin,xMin],[yMin,yMax],'m--');
    hmax=plot([xMax,xMax],[yMin,yMax],'m-.');

    legendStrings={quantityString, 'left bound', 'right bound'};
    legend(legendStrings,'Location','NorthWest');
    handleText=text(mean(x),mean(y),sprintf('rms=%.3f%%',100*rms));
%    axis([-inf,inf,min(floor((y*10))/10),max(ceil((y*10))/10)]);
    axis([valueMin,valueMax,yMin,yMax]);
    
    if interactive==2 || (interactive==4 && rms>rmsCrit)
        disp('INFO (getSlope): indicate lower/upper bound; press <enter> to confirm or <ESC> to run in batch mode.')
        button=1;
        while(button)
            [xValue,yValue,button]=ginput(1);
            hold on
            if ~isempty(button)
                if(button==27) % ESC key
                    interactive=0; % temporary setting to 0
                    disp('INFO (getSlope): running in batch mode')
                    break;
                end

                switch button
                    case 1
                        xMin=xValue;
                        delete(hmin);
                        hmin=plot([xMin,xMin],[yMin,yMax],'m--');
                    case 3
                        xMax=xValue;
                        delete(hmax);
                        hmax=plot([xMax,xMax],[yMin,yMax],'m-.');
                    otherwise
                        disp('WARNING (getSlope): unknown button/key pressed')
                end
                if xMin>xMax
                    xStore=xMax;
                    xMax=xMin;
                    xMin=xStore;
                end
                [coeff,fit,rms,rSquared,error]=fitSlope(type,xMin,xMax,x,y);
                delete(hfit);
                hfit=plot(x,fit,'r-.');
                axis([valueMin,valueMax,yMin,yMax]);
                
                delete(handleText);
                handleText=text(mean(x),mean(y),sprintf('rms=%.3f%%',100*rms));

            end
        end
        if xMin>xMax
            xDummy=xMax;
            xMaxxMin;
            xMin=xDummy;
        end
    end
    close(hfig);
    fitParameters=coeff;
    
end
