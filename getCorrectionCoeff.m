%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Generating signal correction parameters
% based on slope-matching: CO2->O2, O2->MMss
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% 3.0, 26. Juli 2012
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters]=getCorrectionCoeff(withRecalibration,filePath,parameters,type,graphState)

    fprintf('parameter calibration for signal processing with data from %s\n\n',filePath);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set general parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    graphState  =   parameters.Simulation.graphState;
    response    =   parameters.Operation.response;
    dt          =   parameters.Simulation.dt;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reading data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [header, data, signal]  =   readSignalSpirowareMassSpec(filePath,1,parameters);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interactive cropping of data set for offset calibration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if parameters.Simulation.torontoFile
        type = 1;
    end
    
    switch(type)
        case 1
            [time.min,time.max]	=   cropDataSpirowareMassSpec(signal,parameters,'calibration interval');
        case 0
            time	=	parameters.Calibration.time;
        case -1
            time.min    =   min(signal.ts);
            time.max	=   max(signal.ts);
        case 2
            fprintf('\nautomatic calibration\n')
            time.max    =   find(signal.MMss > (max(signal.MMss)+min(signal.MMss) )/2, 1, 'first'); 
            if isempty(time.max)
                h = getOrMakeFigure('Faulty O2 signal');
                set(0, 'currentfigure', h);
                plot(signal.O2);
                error(['No Washin/Washout detected. Aborting processing of File: ' filePath]);
            end
            
            [breathIndexes,breathTimes]=breathDetection(signal.ts(1:time.max),signal.Iv(1:time.max),parameters);            
                                     
            %if at least 7 occurences of an edge are found do the
            %calibration automatically, else do it manually
            if length(breathIndexes)>=3
                time.min = max(breathIndexes(max(end-7,1))-50,1);
                time.max = signal.ts(time.max-10);
                time.min = signal.ts(time.min);
            else
                warning('No 3 breaths found');
               [time.min,time.max] = cropDataSpirowareMassSpec(signal,parameters,'calibration interval');
            end
            
            visualCheck=1;
            if visualCheck && graphState
                plotVisualCheck(signal, time);
            end   
    end

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Crop signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    imin    =   min(find(signal.ts>=time.min));
    imax    =   max(find(signal.ts<=time.max));
    
    fields              =   fieldnames(signal);
    for i=1:length(fields)
        if length(signal.(fields{i}))~=1
            signal.(fields{i})  =   signal.(fields{i})(imin:imax);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setting default values and determination of coarse delays for O2 and MMss
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~response
        tauCO2          =   0.05;                           % start value response accommodation CO2->O2
    else
        tauCO2          =   parameters.Simulation.tauFix;	% fixed response accommodation O2->CO2
    end
    kCO2                =   1.0;
    kO2                 =   1.0;
    kMMss               =   1.0;
    [delayO2,delayMMss] =	getDelayCoarse(signal,parameters);
    
    delayCorrection     =   tauCO2;
    delayO2             =   delayO2-delayCorrection;
    delayMMss           =   delayMMss-delayCorrection;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Work around to treat negative delay issues
    % Basically, the O2 and the MMss signal are retarded by -min(delayO2,delayMMss)*1.5
    % in order to be able to apply the filter optimization mechanism
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minDelay = min(delayO2,delayMMss);
    if minDelay < 0
        parameters.Simulation.nAdvancingSamples =   ceil(-1.5*minDelay/dt);
    else
        parameters.Simulation.nAdvancingSamples =   0;  % no artifical advancing necessary
    end
    
    signalTemp	=   signalAdvance(signal,parameters);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Optimize filter coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coefficients        =   [delayO2,delayMMss,tauCO2,kCO2,kO2,kMMss];
    coefficientScales   =   [1,      1,        0.01,  1,   1,  1    ];
    
    %optimizingSet       =   [1,2,3,4,6];                % coefficient indices to optimize (MMss)
    %optimizingSet       =   [1,2,3,4];                  % coefficient indices to optimize (N2)
    if ~response
        optimizingSet	=   [1,2,3,4,6];                % coefficient indices to optimize (MMss+N2, with tauCO2)
    else
        optimizingSet	=   [1,2,4,6];                  % coefficient indices to optimize (MMss+N2, without tauCO2)
    end
    
    coefficientsOptimal	=   optimizeFilters(optimizingSet,coefficientScales,coefficients,signalTemp,parameters);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Enforce absence of recalibration
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~withRecalibration
        coefficientsOptimal([4,5,6])	=   1;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Returning calibration values in parameters data structure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parameters.Calibration.coefficientsOptimal  =   coefficientsOptimal;
    parameters.Calibration.time                 =   time;
    parameters.Calibration.filePath             =   filePath;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Apply filters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fprintf('\nAppliying signal correction to source file\n');
    signal  =   signalCorrection(signal,parameters,graphState);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Self consistency test information                                                         
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    meanN2  =   mean(signal.N2Poi3);
    stdN2   =   std(signal.N2Poi3);
    fprintf('meanN2=%g, stdN2=%g, errorN2=%g %%\n' ,meanN2,stdN2,stdN2/meanN2*100);
end

function plotVisualCheck(signal, time)
    
    h = getOrMakeFigure('Visual Calibration Check');
    set(0, 'currentfigure', h);
    plot(signal.ts, signal.CO2, signal.ts, signal.O2, signal.ts, signal.MMss);
    legend('CO2', 'O2', 'MMss');
    hold on
    line([time.min, time.min], [min(signal.O2), max(signal.O2)], 'LineWidth',2, 'Color', 'black');
    line([time.max, time.max], [min(signal.O2), max(signal.O2)], 'LineWidth',2, 'Color', 'black');
    line([time.min, time.max], [max(signal.O2)*0.9, max(signal.O2)*0.9], 'LineWidth', 1, 'Color', 'black', 'LineStyle','--');
    hold off
    
end
