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
function [signal,breathIndex]=signalCorrection(signal,parameters,flag)

    torontoFile     =   parameters.Simulation.torontoFile;
    bFile           =   parameters.Simulation.bFile;
    nOrder          =   parameters.Operation.filterOrder;
    filterCutoff    =   parameters.Operation.filterCutoff;
    dt              =   parameters.Simulation.dt;
    graphState      =   parameters.Simulation.graphState;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % get optimal correction coefficients
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    coefficients    =   parameters.Calibration.coefficientsOptimal; % optimal correction coeffients
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % apply signal correction filter and calculate MMcalc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if torontoFile 
        signal.IvFilter         =	signal.Iv;          % dummy filtering (for the time beeing)
        signal.delayLineFilter	=   signal.delayLine;	% dummy filtering
        signal.delayTempFilter	=   signal.delayTemp;	% dummy filtering
    elseif bFile
        signal.CO2Filter   	=   signal.CO2;             % dummy filtering (for the time beeing)
        signal.O2Filter   	=   signal.O2;              % dummy filtering (for the time beeing)
        signal.MMssFilter   =   signal.MMss;            % dummy filtering (for the time beeing)
        signal.MMssDelay    =   signal.MMss;            % dummy filtering (for the time beeing)
        signal.IvFilter     =   signal.Iv;              % dummy filtering (for the time beeing)
        signal.IvDelay      =   signal.Iv;              % dummy filtering (for the time beeing)
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % manually  change correction coefficients
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        deltaCO2O2              =   parameters.Operation.deltaCO2O2;
        deltaCO2MMss            =   parameters.Operation.deltaCO2MMss;
        dTauCO2                 =   parameters.Operation.dTauCO2;
        if isempty(coefficients)
           coefficients = zeros(1, 6); 
        end
        coefficients            =   coefficients + [ deltaCO2O2 deltaCO2MMss dTauCO2 0 0 0 ];
                
        signal                  =   filterDelayTemp(signal,parameters);
        nOrder                  =   1;
        signal                  =   filterRawSignal(signal, nOrder, filterCutoff, dt); 
        signal                  =   filterDelayIV(signal,parameters);
        signal                  =   dataDelay(signal,parameters);  
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate the gas concentrations in the signal
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signal  =   signalCalc(signal,parameters);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Breath segmentation and calculation of volume (breath detected variant)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [breathIndex,breathTimes]=breathDetection(signal.ts,signal.IvEffBTPS,parameters);
    signal.breathIndex = breathIndex;
    signal.breathTimes = breathTimes;
    
    VolEffBTPS      =   cumtrapz(signal.ts,signal.IvEffBTPS);	% calculating the "integrated" volumes
    VolBTEffBTPS	=   zeros(size(VolEffBTPS));
    breathTimes = zeros(size(length(breathIndex)-1));
    for i=1:length(breathIndex)-1
        breathTimes(i)=signal.ts(breathIndex(i+1));
        gapVol=VolEffBTPS(breathIndex(i));
        VolBTEffBTPS(breathIndex(i):breathIndex(i+1)-1)=VolEffBTPS(breathIndex(i):breathIndex(i+1)-1)-gapVol;
    end
    signal.VolEffBTPS	=   VolEffBTPS;
    signal.VolBTEffBTPS	=   VolBTEffBTPS;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % plain integration of volume
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signal.VolFilter	=   cumtrapz(signal.ts,signal.IvDelay+signal.Ivss);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting results
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag~=0 && ~bFile && graphState
        plotCorrectedSignal( signal, parameters, breathTimes, flag );
    end
end