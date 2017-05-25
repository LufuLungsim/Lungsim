%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Interactive calculation of slope, works also for batch operation
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 12. Dezember 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [capnoIndices,interactive,errorInfo]=getCapno(minX,maxX,minY,maxY,x,y,flow,interactive,rmsCrit)

    capnoIndices.isValid	= 1;
    errorInfo=0;
    y=y-y(1);

    xMin=min(x);
    yMin=min(y);
    xMax=max(x);
    yMax=max(y);
    
    xMinIII=minX*xMax;
    xMaxIII=maxX*xMax;
    indicesIII=find(x>xMinIII & x<xMaxIII);
    yMinIII=min(y(indicesIII));
    yMaxIII=max(y(indicesIII));

    yMinII =minY*yMax;
    yMaxII =maxY*yMax;
    indicesII =find(y>yMinII  & y<yMaxII );
    xMinII =min(x(indicesII));
    xMaxII =max(x(indicesII));
    
    hfig=figure(902);
    if interactive>1    % adjusting figure size if interactive
        figureSize=get(hfig,'Position');
        height=figureSize(4);
        figureSize(4)=1.6*height;
        figureSize(2)=figureSize(2)-0.6*height;
        set(hfig,'Position',figureSize);
    end
    set(hfig,'Name',strcat('Capno index determination'));
    set(hfig,'NumberTitle','off');
    
    subplot(3,1,2)
    plot(x,flow,'b');
    title('volume flow plot');
    xlabel('Vol / m^3');
    ylabel('Iv / m^3/');
    axis([xMin,xMax,-inf,inf]);
    
    subplot(3,1,1)
    hold off
    plot(x,y,'b');
    title({'Indicate lower/upper bound, <x>,<d> dismiss','','<enter> accepts breath, <ESC> runs batch mode'});
    xlabel('Vol / m^3');
    ylabel('CO2');
    hold on
    
    if xMaxII>xMinIII
        capnoIndices.isValid   = 0;
        capnoIndices.fCO2E     = 0;
        capnoIndices.CO2et     = 0;
        capnoIndices.coeffII   = 0;
        capnoIndices.coeffIII  = 0;
        capnoIndices.fowlerDead= 0;
        if interactive>1
            h = warndlg('Breath cannot be analysed (wrong II/III-intervals)','getCapno');
            uiwait(h);
        else
            fprintf('Warning (getCapno): Breath cannot be analysed (wrong II/III-intervals)\n')
        end
        close(hfig);
        return;
    end
    
    yDiff=yMax-yMin;
    yMin=yMin-yDiff/10;
    yMax=yMax+yDiff/10;
    
    fowlerLeft=cumtrapz(x,y);
    vCO2E=fowlerLeft(end);

    [coeffIII,fitIII,rmsIII,rSquaredIII,errorInfo]=fitSlope('lin',xMinIII,xMaxIII,x,y);
    if errorInfo
        capnoIndices.isValid   = 0;
        capnoIndices.fCO2E     = 0;
        capnoIndices.CO2et     = 0;
        capnoIndices.coeffII   = 0;
        capnoIndices.coeffIII  = 0;
        capnoIndices.fowlerDead= 0;
        if interactive>1
            h = warndlg('Breath cannot be analysed (slopeIII)','getCapno');
            uiwait(h);
        else
            fprintf('Warning (getCapno): Breath cannot be analysed (slopeIII)\n')
        end
        close(hfig);
        return;
    end
    hfitIII=plot(x,fitIII,'r-','LineWidth',2.0);
    hminIII=plot([xMinIII,xMinIII],[yMin,yMax],'r-.');
    hmaxIII=plot([xMaxIII,xMaxIII],[yMin,yMax],'r-.');
    handleTextIII=text((xMinIII+xMaxIII)/2,(yMinII+yMaxII)/4*3,sprintf('sIII=%.3f %%/l\nsnIII=%.3f 1/l\nrms=%.3f %%',coeffIII(2)*1e2/1e3,coeffIII(2)/(vCO2E/xMax)/1e3,rmsIII*1e2),'Color','r');

    [coeffII ,fitII ,rmsII ,rSquaredII ,errorInfo]=fitSlope('lin',xMinII,xMaxII,x,y);
    if errorInfo
        capnoIndices.isValid   = 0;
        capnoIndices.fCO2E     = 0;
        capnoIndices.CO2et     = 0;
        capnoIndices.coeffII   = 0;
        capnoIndices.coeffIII  = 0;
        capnoIndices.fowlerDead= 0;
        if interactive>1
            h = warndlg('Breath cannot be analysed (slopeII)','getCapno');
            uiwait(h);
        else
            fprintf('Warning (getCapno): Breath cannot be analysed (slopeII)\n')
        end
        close(hfig);
        return;
    end
    hfitII =plot(x,fitII,'m-','LineWidth',2.0);
    hminII =plot([xMinII ,xMinII ],[yMin,yMax],'m-.');
    hmaxII =plot([xMaxII ,xMaxII ],[yMin,yMax],'m-.');
    handleTextII =text((xMinII +xMaxII )/2,(yMinII+yMaxII)/4*3,sprintf('sII=%.3f %%/l\nsnII=%.3f 1/l\nrms=%.3f %%',coeffII(2)*1e2/1e3,coeffII(2)/(vCO2E/xMax)/1e3,rmsII*1e2),'Color','m');

    axis([xMin,xMax,0,yMax]);
    
    fowlerRight=cumtrapz(x,y-fitIII);
    fowlerRight=fowlerRight-fowlerRight(end);
    fowlerDead=fowlerLeft-fowlerRight;
    indexDead=find(fowlerDead<0);
    if isempty(indexDead)
        errorInfo=1;
        capnoIndices.isValid   = 0;
        capnoIndices.fCO2E     = 0;
        capnoIndices.CO2et     = 0;
        capnoIndices.coeffII   = 0;
        capnoIndices.coeffIII  = 0;
        capnoIndices.fowlerDead= 0;
        if interactive>1
            h = warndlg('Breath cannot be analysed (dead space)','getCapno');
            uiwait(h);
        else
            fprintf('Warning (getCapno): Breath cannot be analysed (dead space)\n')
        end
        close(hfig);
        return;
    end
    xFowlerDead=fowlerDead(indexDead(end)+1)/(fowlerDead(indexDead(end)+1)-fowlerDead(indexDead(end)))*x(indexDead(end))-fowlerDead(indexDead(end))/(fowlerDead(indexDead(end)+1)-fowlerDead(indexDead(end)))*x(indexDead(end)+1);
    yFowler=y(indexDead(end));
    
    hDeadII =plot([xFowlerDead ,xFowlerDead ],[yMin,yFowler],'m-');
    hDeadTextII =text(xFowlerDead*0.8,yFowler/5,sprintf('Vdead=%.1f ml',1e6*xFowlerDead ));
    
    subplot(3,1,3)
    hold off
    plot(x,zeros(size(x)),'k--');
    hold on
    hFowlerLeft =plot(x,fowlerLeft, 'b');
    hFowlerRight=plot(x,fowlerRight,'r');
    hFowlerDead =plot(x,fowlerDead, 'g');
    xlabel('Vol / m^3');
    ylabel('Fowler integrals');
    axis([xMin,xMax,-vCO2E/20,vCO2E/20]);
    
    if interactive==2 || (interactive==4 && (rmsII>rmsCrit || rmsIII>rmsCrit ))
        disp('INFO (getCapno): indicate lower/upper bound.\nPress <x> or <d> to dismiss, <enter> to accept or <ESC> to run in batch mode.')
        button=1;
        while(button)
            [xValue,yValue,button]=ginput(1);
            hold on
            if ~isempty(button)
                if(button==27) % ESC key
                    interactive=0; % temporary setting to 0
                    disp('INFO (getCapno): running in batch mode')
                    break;
                elseif(button==100 || button==120)
                    disp('INFO (getCapno): dismiss this breath')
                    capnoIndices.isValid   = 0;
                    break;
                end

                subplot(3,1,1)
                switch button
                    case 1
                        if abs(xValue-xMinII) < abs(xValue-xMinIII)
                            xMinII=xValue;
                            delete(hminII);
                            hminII =plot([xMinII ,xMinII ],[yMin,yMax],'m-.');
                        else
                            xMinIII=xValue;
                            delete(hminIII);
                            hminIII =plot([xMinIII ,xMinIII ],[yMin,yMax],'r-.');
                        end
                    case 3
                        if abs(xValue-xMaxII) < abs(xValue-xMaxIII)
                            xMaxII=xValue;
                            delete(hmaxII);
                            hmaxII =plot([xMaxII ,xMaxII ],[yMin,yMax],'m-.');
                        else
                            xMaxIII=xValue;
                            delete(hmaxIII);
                            hmaxIII =plot([xMaxIII ,xMaxIII ],[yMin,yMax],'r-.');
                        end
                    otherwise
                        disp('WARNING (getCapno): unknown button/key pressed')
                end
                if xMinII>xMaxII
                    xStore=xMaxII;
                    xMaxII=xMinII;
                    xMinII=xStore;
                end
                if xMinIII>xMaxIII
                    xStore=xMaxIII;
                    xMaxIII=xMinIII;
                    xMinIII=xStore;
                end
                [coeffII,fitII,rmsII,rSquaredII ,errorInfo]=fitSlope('lin',xMinII,xMaxII,x,y);
                delete(hfitII);
                hfitII =plot(x,fitII,'m-.','LineWidth',2.0);
                delete(handleTextII);
                handleTextII =text((xMinII +xMaxII )/2,(yMinII+yMaxII)/4*3,sprintf('sII=%.3f %%/l\nsnII=%.3f 1/l\nrms=%.3f %%',coeffII(2)*1e2/1e3,coeffII(2)/(vCO2E/xMax)/1e3,rmsII*1e2),'Color','m');
                
                [coeffIII,fitIII,rmsIII,rSquaredIII,errorInfo]=fitSlope('lin',xMinIII,xMaxIII,x,y);
                delete(hfitIII);
                hfitIII=plot(x,fitIII,'r-.','LineWidth',2.0);
                delete(handleTextIII);
                handleTextIII=text((xMinIII+xMaxIII)/2,(yMinII+yMaxII)/4*3,sprintf('sIII=%.3f %%/l\nsnIII=%.3f 1/l\nrms=%.3f %%',coeffIII(2)*1e2/1e3,coeffIII(2)/(vCO2E/xMax)/1e3,rmsIII*1e2),'Color','r');
                
                axis([xMin,xMax,0,yMax]);
                
                fowlerLeft=cumtrapz(x,y);
                vCO2E=fowlerLeft(end);
                fowlerRight=cumtrapz(x,y-fitIII);
                fowlerRight=fowlerRight-fowlerRight(end);
                fowlerDead=fowlerLeft-fowlerRight;
                indexDead=find(fowlerDead<0);
                xFowlerDead=fowlerDead(indexDead(end)+1)/(fowlerDead(indexDead(end)+1)-fowlerDead(indexDead(end)))*x(indexDead(end))-fowlerDead(indexDead(end))/(fowlerDead(indexDead(end)+1)-fowlerDead(indexDead(end)))*x(indexDead(end)+1);
                yFowler=y(indexDead(end));
                
                delete(hDeadII);
                delete(hDeadTextII);
                hDeadII =plot([xFowlerDead ,xFowlerDead ],[yMin,yFowler],'m-');
                hDeadTextII =text(xFowlerDead*0.8,yFowler/5,sprintf('Vdead=%.1f ml',1e6*xFowlerDead ));
                
                subplot(3,1,3)
                delete(hFowlerLeft);
                delete(hFowlerRight);
                delete(hFowlerDead);
                hFowlerLeft =plot(x,fowlerLeft,'b');
                hFowlerRight=plot(x,fowlerRight,'r');
                hFowlerDead =plot(x,fowlerDead,'g');                
            end
        end
        if xMinII>xMaxII
            xStore=xMaxII;
            xMaxII=xMinII;
            xMinII=xStore;
        end
        if xMinIII>xMaxIII
            xStore=xMaxIII;
            xMaxIII=xMinIII;
            xMinIII=xStore;
        end
    end
    close(hfig);
    
    capnoIndices.fCO2E     = vCO2E/xMax;
    capnoIndices.CO2et     = median(y(round(0.98*length(y)):end));
    capnoIndices.coeffII   = coeffII(2)/capnoIndices.fCO2E;
    capnoIndices.coeffIII  = coeffIII(2)/capnoIndices.fCO2E;
    capnoIndices.fowlerDead= xFowlerDead;

end
