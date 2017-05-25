%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to filter signals
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.4, 20. Oktober 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal]=signalFilter(signal,parameters,coefficients)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % fixed values of O2, CO2 and MMss
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    dt                      =	parameters.Simulation.dt;                   % time step
    dtAdvance               =	parameters.Simulation.nAdvancingSamples*dt; % optional advancing time
    Ivss0                   =	parameters.Operation.sampleFlow;            % nominal sample flow
    lag                     =	parameters.Operation.maximalLag;            % maximal lag of signals
    verb                    =	parameters.Simulation.verb;                 % verbosity level
    slow                    =	parameters.Simulation.slowState;            % slowing down response
    aFile                   =	parameters.Simulation.aFile;                % file type chosen
    delayFlow               =	parameters.Calibration.delayFlow;           % non-detected delay of flow and species signals
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filters for O2 related corrections
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    delayO2                 =   coefficients(1);                            % constant delay ms->ss
    delayMMss               =   coefficients(2);                            % constant delay ms->ss
    tauCO2                  =   coefficients(3);                            % response accommodation CO2->O2
    kCO2                    =   coefficients(4);                            % sensitivity CO2 sensor
    kO2                     =   coefficients(5);                            % sensitivity O2  sensor
    kMMss                   =   coefficients(6);                            % sensitivity MMss sensor (temperature coefficient?)
    
    maxDelay                =   max(delayO2,delayMMss);                     % maximal delay
    
    varyingMaxDelay         =   (dtAdvance+maxDelay)*Ivss0./abs(signal.IvssSmooth);
    varyingO2Delay          =   (dtAdvance+maxDelay-delayO2)*Ivss0./abs(signal.IvssSmooth);
    varyingMMssDelay        =   (dtAdvance+maxDelay-delayMMss)*Ivss0./abs(signal.IvssSmooth);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delays (translations in time)
    % new setup: CO2 in ms, approximately synchronus with flow signal. The flow signal
    %            is therefore treated according to the CO2 signal (+delayFlow).
    % old setup: CO2 in ss, but it is assumed that by exporting from wbreath, the delays
    %            of the species signals relative to the flow have already been corrected.
    %            Only interspecies corrections are performed
    % Conclusion: old and new setup are treated in the same way w.r.t. constant
    % delay but for different reasons, namely old: delay is already corrected at
    % export, new: flow is treated synchronously to CO2
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    bFlow                   =	filterLagrange((varyingMaxDelay+delayFlow)/dt,3);	% constant delay Flow
    bCO2                    =   filterLagrange(varyingMaxDelay/dt,3);            	% constant delay CO2
    bO2                     =   filterLagrange(varyingO2Delay/dt,3);             	% constant delay O2
    bMMss                   =   filterLagrange(varyingMMssDelay/dt,3);           	% constant delay MMss
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Delay corrections
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signal.CO2Filter        =   filterFIR(bCO2,signal.CO2);                         % delay CO2 signal
    signal.IvFilter         =   filterFIR(bFlow,signal.Iv);                         % delay on Iv signal
    signal.IvssFilter       =   signal.Ivss;
    signal.IvssSmoothFilter =   signal.IvssSmooth;
    signal.IvssSmooth1Filter=   signal.IvssSmooth1;
    signal.delayLineFilter	=   filterFIR(bFlow,signal.delayLine);                  % delay on delayLine signal
    signal.delayTempFilter	=   filterFIR(bFlow,signal.delayTemp);                  % delay on delayTemp signal
    
    signal.MMmsFilter       =   filterFIR(bCO2,signal.MMms);                        % delay MM signal
    signal.O2Filter         =   filterFIR(bO2,signal.O2);                           % delay O2 signal
    signal.MMssFilter       =   filterFIR(bMMss,signal.MMss);                       % delay MMss signal
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Sensor response corrections: first order IIR filter to accommodate CO2/Iv to
    % O2/MMss signal
    % There are to variants depending on the setup chosen:
    %   a) new setup: CO2-sensor is the fastest species sensor (flow sensors iseven faster)
    %   b) old setup: CO2-sensor is slowest sensor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if aFile
        if ~slow
            fButterworth	=   8.0;  % was 8.0                                 	% Butterworth critical frequency
            tauCO2          =   tauCO2+1/(fButterworth*2*pi);                      	% tauCO2 modification to accomodate slowing of filter!
        end
        pCO2                =   1/tauCO2;
        dCO2                =   [1,pCO2];                                           % denominator of 1st order continuous filter
        nCO2                =   dCO2(end);                                          % nominator of second order continuous filter
        [denCO2,nomCO2]    	=   filterBilinear(nCO2,dCO2,0,dt);                     % second order discret filter by bilinear transformation
        
        if ~slow % accelerating resonse of O2, MMss
            [a,b]               =   filterButterworth(fButterworth,1,dt);           % discret butterworth filter by bilinear transformation
            signal.O2Butter     =   filterIIR(b,a,signal.O2Filter,[]);              % filtering signals to prevent excessive noise
            signal.MMssButter   =   filterIIR(b,a,signal.MMssFilter,[]);
            
            signal.O2Filter1	=   filterIIR(denCO2,nomCO2,signal.O2Butter,[]);    % inverse response filter
            signal.MMssFilter1	=   filterIIR(denCO2,nomCO2,signal.MMssButter,[]);
            
            if verb
                figure(701)
                subplot(2,1,1)
                plot(signal.ts,[signal.O2Filter,  signal.O2Butter,  signal.O2Filter1  ]);
                subplot(2,1,2)
                plot(signal.ts,[signal.MMssFilter,signal.MMssButter,signal.MMssFilter1]);
            end
            
            signal.O2Filter     =   signal.O2Filter1;                               % mapping to correct signals
            signal.MMssFilter   =   signal.MMssFilter1;
            
        else    % slowing down response of CO2, flow
            signal.CO2Filter        =   filterIIR(nomCO2,denCO2,signal.CO2Filter,[]);       % response filter
            signal.IvFilter         =   filterIIR(nomCO2,denCO2,signal.IvFilter,[]);       	% response filter applied to Iv signal (for synchronization)
            signal.delayLineFilter	=   filterIIR(nomCO2,denCO2,signal.delayLineFilter,[]);	% response filter applied to delayLine signal (for synchronization)
        end
    else    % old setup
        if ~slow
            fButterworth	=   8.0;  % was 8.0                                     % Butterworth critical frequency
            tauCO2          =   tauCO2+1/(fButterworth*2*pi);                       % tauCO2 modification to accomodate slowing of filter!
        end
        pCO2                =   1/tauCO2;
        dCO2                =   [1,pCO2];                                           % denominator of 1st order continuous filter
        nCO2                =   dCO2(end);                                          % nominator of second order continuous filter
        [denCO2,nomCO2]    	=   filterBilinear(nCO2,dCO2,0,dt);                     % second order discret filter by bilinear transformation
        
        if ~slow	% accelerating CO2 response
            [a,b]               =   filterButterworth(fButterworth,1,dt);           % discret butterworth filter by bilinear transformation
            signal.CO2Butter    =   filterIIR(b,a,signal.CO2Filter,[]);             % filtering signals to prevent excessive noise
            
            signal.CO2Filter1	=   filterIIR(denCO2,nomCO2,signal.CO2Butter,[]);   % inverse response filter
            
            if verb
                figure(702)
                plot(signal.ts,[signal.CO2Filter,  signal.CO2Butter,  signal.CO2Filter1  ]);
            end
            
            signal.CO2Filter    =   signal.CO2Filter1;                              % mapping to correct signals
            
        else        % slowing down response of O2, MMss, Iv, etc.
            signal.O2Filter         =   filterIIR(nomCO2,denCO2,signal.CO2Filter,[]);       % response filter
            signal.MMssFilter       =   filterIIR(nomCO2,denCO2,signal.MMssFilter,[]);      % response filter
            signal.IvFilter         =   filterIIR(nomCO2,denCO2,signal.IvFilter,[]);       	% response filter applied to Iv signal (for synchronization)
            signal.delayLineFilter	=   filterIIR(nomCO2,denCO2,signal.delayLineFilter,[]);	% response filter applied to delayLine signal (for synchronization)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Gain correction for CO2, O2 and MMss
    % O2 gain: could be implemented to fit kO2 to given O2 plateau!!
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    signal.CO2Filter   	=   (signal.CO2Filter -0)*kCO2+0;	% CO2 response filter (gain only)
    signal.O2Filter   	=   (signal.O2Filter  -0)*kO2+0;	% O2 response filter (gain only)
    signal.MMssFilter   =   (signal.MMssFilter-0)*kMMss+0;	% MMss response filter (gain only)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Signal truncatoin (cut away inconsistent signal portion)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nTruncation         =   floor(lag/dt);
    fields              =   fieldnames(signal);
    for i=1:length(fields)
        signal.(fields{i})  =   signal.(fields{i})(nTruncation+1:end);
    end
end