 %
% Manual Outlier Detection for LCI data
% Lets the user select datapoints to be integrated in the calculations of
% the LCI by FIT
%
function validData = manualOutlierDetection(x,y,validData,signalCorrected,gasByName,parameters)

    k=0;    
    hfig=figure(3562);
    set(hfig,'Name',sprintf('Manual data validation for %s',gasByName));    
    subplotA = drawFigureDescription(hfig, x,y, signalCorrected, validData, parameters);
    set(hfig,'Visible','on');
    subplot(3,1,1);
    
    figureSize=get(hfig,'Position');    % rescaling graph if necessary (for better visibility)
    height=figureSize(4);
    width=figureSize(3);
    if width<1000
        figureSize(4)=2.1*height;
        figureSize(2)=figureSize(2)-1.1*height;
        figureSize(3)=2.1*width;
        figureSize(1)=figureSize(1)-1.1*width;
        set(hfig,'Position',figureSize);
    end

    
    %k is set by waitforbuttonpress which returns 1 if a key is pressed and
    %2 if a mouse is clicked, so do while k is not set by the program
    while k ~= 3
        k=waitforbuttonpress;
        if k == 1                   %if keypress k=1
            if strcmp(get(gcf,'currentcharacter'),'v');
                validData=updateLCIFit(hfig, subplotA, x,y,validData, signalCorrected);                
            elseif strcmp(get(gcf,'currentcharacter'),'e'); 
                k=3;
            end           
        end
   end  
   %close(hfig);    
end


%
% select datapoints to exclude or include for calculation of the LCI by fit function
%
function validData=updateLCIFit(hfig, subplotA, x,y,validData, signal)
    set(hfig, 'currentaxes', subplotA);
    [xm1 ym1] = ginput(1);  %Get first corner of the selection rectangle
    hold on;
    line([xm1 xm1],         [min(y)*0.9 max(y)*1.1]); %Draw the borders of the selection rectangles
    line([0  max(x)*1.1],   [ym1 ym1]);
    hold off;
    [xm2 ym2] = ginput(1);

    %flip the selectedPoints in the selected area to include/exclude the datapoint                          
    selectedPoints=x>min([xm1 xm2])&x<max([xm1 xm2])&y>min([ym1 ym2])&y<max([ym1 ym2]);
    validData(selectedPoints)=1-validData(selectedPoints); 
    
    [params,rms,r2,handle] = multiExponentialFit(x(validData==1),y(validData==1));

    plotLCIFunction(hfig, subplotA, x, y, validData, handle);
end

%
% plot the LCI by fit function and the corresponding dataPoints
%
function plotLCIFunction(hfig, axes, x, y, validData, handle)
    
    set(hfig, 'currentaxes', axes);
    semilogy(x,handle(x), '--r');                                   %Plot fit function on relevant data
    hold(axes, 'on');
    semilogy(x(validData==1), y(validData==1), '.b', 'MarkerSize', 20);               %Plot datapoints used to find the fit function
    semilogy(x(validData==0),y(validData==0), 'ob');
    line([0  max(x)*1.1], [0.025,0.025], 'LineStyle', '--', 'Color', 'g');
    axis([0  max(x)*1.1 min(y)*0.9 max(y)*1.1]);
    
    text(x,y,strcat(' ', num2str([1:length(x)]')));
    title({'Manually validate the data to be analized. '...
           'Filled circles are included, empty circles are ignored.'... 
           'Press "v" to select data to be in/excluded'...
           'or "e" to continue'});
    xlabel('CEV/l');
    ylabel('Cet(norm)/%');
    axis([0  max(x)*1.1 min(y)*0.9 max(y)*1.1]);
    hold(axes, 'off');
end

%
% Plot initial Figure containing current LCI by Fit function, gasvolume
% plot and gasflow plot
%
function subplotA = drawFigureDescription(hfig, x,y, signal, validData, parameters)
    torontoFile=parameters.Simulation.torontoFile;
    
    [MBWPhase, adapdedXAxis, breathIndices] = getMBWPhaseAndAdapdedAxis(signal, x, parameters);

    subplotA=subplot(3,1,1);
    [params,rms,r2,handle] = multiExponentialFit(x(validData==1),y(validData==1));
    plotLCIFunction(hfig, subplotA, x, y, validData, handle)
    
    if torontoFile
        subplotB=subplot(3,1,2);
        set(hfig, 'currentaxes', subplotB);
        plot(adapdedXAxis,signal.CO2(MBWPhase)*100, adapdedXAxis,signal.He(MBWPhase)*100, adapdedXAxis,signal.SF6(MBWPhase)*100);
        title('molar fractions at POI');
        xlabel('Adapded timeaxis to top plot');
        ylabel('x / %');
        legend('CO2', 'He', 'SF6');
        axis([-inf,inf,-inf,inf]);
        annotatePlot(hfig, subplotB, breathIndices, x, 0, 10, 70, signal.VolBTEffBTPS(MBWPhase)*1e6);
    else
        subplotB=subplot(3,1,2);
        set(hfig, 'currentaxes', subplotB);
        plot(adapdedXAxis,signal.CO2Poi1(MBWPhase)*100, adapdedXAxis,signal.O2Poi1(MBWPhase)*100, adapdedXAxis,signal.N2Poi1(MBWPhase)*100,...
            adapdedXAxis,signal.ArPoi1(MBWPhase)*100);
        title('molar fractions at POI1');
        xlabel('Adapded timeaxis to top plot');
        ylabel('x / %');
        legend('CO2', 'O2', 'N2', 'Ar');
        axis([-inf,inf,-inf,inf]);
        annotatePlot(hfig, subplotB, breathIndices, x, 0, 100, 90, signal.VolBTEffBTPS(MBWPhase)*1e6);
    end
    
    subplotC=subplot(3,1,3);
    set(hfig, 'currentaxes', subplotC);
    plot(adapdedXAxis,signal.IvEffBTPS(MBWPhase)*1e6,adapdedXAxis,signal.VolBTEffBTPS(MBWPhase)*1e6);
    title('flow related data');
    xlabel('Adapded timeaxis to top plot');
    ylabel('I_V_,_B_T_P_S, V_B_T_P_S / ml/s, ml');
    legend('IvEffBTPS', 'VolBTEffBTPS');
    axis([-inf,inf,-inf,inf]);
    annotatePlot(hfig, subplotC, breathIndices, x, -500, 500, 10, signal.VolBTEffBTPS(MBWPhase)*1e6);
    linkaxes([subplotA, subplotB, subplotC],'x');
end

%
% Annotate breaths
%
function annotatePlot(hfig, axes, breathIndices, x, xmin, xmax, annotationHight, signal)
    set(hfig, 'currentaxes', axes);
    hold on;
    indexFirstBreath = find(signal==0, 1);
    increase=signal(indexFirstBreath)-signal(indexFirstBreath+10);
    if increase>0   %Check if in or expirations
        increase=0;
    else
        increase=1;
    end
    
    for time=breathIndices
        plot([time,time],[xmin, xmax],'k-.');
    end  
    inspirations=breathIndices((2-increase):2:end);
    text(inspirations, annotationHight*ones(size(inspirations)), num2str([1:length(inspirations)]'));
    
    hold off;
end

function [MBWPhase, adapdedXAxis, breathIndices] = getMBWPhaseAndAdapdedAxis(signal, x, parameters)
    torontoFile=parameters.Simulation.torontoFile;
    if torontoFile
        %Find signal range where He is below half of washin concentration% (wahshout of He/SF6 has started)
        normConcentration=parameters.Operation.toronto(3);  % assuming equal start of He and SF6
        MBWPhase=find(signal.He<normConcentration/2):length(signal.He);
    else
        %Find signal range where O2 is over 30% (wahshout of N2 has started)
        MBWPhase=find(signal.O2Poi1>0.3,1):length(signal.O2Poi1);
    end
    %Definde adapded x axis
    adapdedXAxis=(1:length(MBWPhase))/length(MBWPhase)*max(x); 
    %Where the signal is reset to 0, a new breath begins
    breathIndices=find(signal.VolBTEffBTPS(MBWPhase(1):signal.breathIndex(end))==0);
    
    %To filter wrongly detected breaths get the distance between the
    %indexes and filter if the difference is below 3
    indexDiff=diff(breathIndices);
    buffer=breathIndices(2:end);
    breathIndices=[breathIndices(1); buffer(indexDiff>10)];
    
    %Scale to alternate timescale
    breathIndices=breathIndices/length(MBWPhase)*max(x);
end
