%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to filter signals to compensate for variable time delay
% this form of the code is in loop form, i.e., not as efficient as MATLAB
% code but easily convertible to other programmeing languages.
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.1, 27. Januar 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal]=dataDelay(signal,parameters)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameters
% Remark: All parameters regarding delays (form device settings in
% Sprioware are put into setParameters. Variables relating to delays with
% respect to the flow start with delay... 
% ...Static refers to constant delay, while ...Dynamic refers to flow
% dependent delays
% Variables starting with delta... are related to pure side stream flow
% delays.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dt              =	parameters.Simulation.dt;                               % time step
delayCO2Static  =   parameters.Calibration.delayCO2;                      	% delay of flow to CO2 signal
delayO2Static   =   parameters.Calibration.delayO2;                      	% delay of flow to O2 signal
delayMMssStatic =   parameters.Calibration.delayMMss;                      	% delay of flow to MMss signal

Ivss0           =   parameters.Operation.sampleFlow;                        % nominal sample flow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% delay evaluation initialization
% Remark: it is assumed that the MMss delay is larger than the O2 delay as
% the MMss sensor is after the O2 sensor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nCompartments	=   256;                                                    % number of compartements used to evaluate delay

deltaO2Static   =   delayO2Static-delayCO2Static;                          	% O2 side stream delay (correct order of magnitude)
deltaMMssStatic =   delayMMssStatic-delayCO2Static;                         % MMss side stream delay (correct order of magnitude)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The setting of the delay line volumina (partVolumeO2 and lineVolume)
% provide only the correct delays, if the value of Ivss0 corresponds to the
% average of the side stream flow during a calibration!
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
partVolumeO2    =	Ivss0*deltaO2Static;                                    % O2 part side stream delay volume
lineVolume      =	Ivss0*deltaMMssStatic;                                  % MMss (total) side stream delay volume

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remark: it is assumed that the MMss delay is larger than the O2 delay as
% the MMss sensor is after the O2 sensor (i.e., lineVolume > partVolumeO2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dVolume         =   lineVolume/nCompartments;                            	% volume increment
leftIndexO2     =   floor(partVolumeO2/dVolume);                            % left index (compartment with volume<partVolumeO2
partIndexO2     =   partVolumeO2/dVolume-leftIndexO2;                       % partial index (non-integer part)
indexO2         =   [leftIndexO2,leftIndexO2+1];                            % index pair containing partVolumeO2
weightsO2       =   [1-partIndexO2,partIndexO2];                            % weights to interpolate to partVolumeO2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% loop over time steps
% This second variant uses fixed time step delaying, where the actual delay
% is rounded to the next integer time step value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
signal.delayMMss        =   zeros(size(signal.ts));                        	% pre declaration of MMss delay signal array
signal.delayO2          =   zeros(size(signal.ts));                        	% pre declaration of O2 delay signal array
signal.CO2Delay         =   zeros(size(signal.ts));
signal.O2Delay          =   zeros(size(signal.ts));
signal.MMssDelay        =   zeros(size(signal.ts));
signal.IvDelay          =   zeros(size(signal.ts));

delayLine               =	ones(1,nCompartments+1)*signal.ts(1);         	% DOFs of the delay line
delayLine               =	delayLine-[0:nCompartments]/(-signal.Ivss(1)/dVolume);	% initialization

nMax                    =   length(signal.ts);

for i=1:nMax
    
    delayLine(1)        =   signal.ts(i);                                   % enforcement of time variing boundary condition
    delayLine           =   singleStepDelay(delayLine,-signal.Ivss(i)*dt/dVolume,nCompartments+1);	% advancing delay line
    deltaMMssDynamic    =   signal.ts(i)-delayLine(end);                    % dynamic MMss delay in side stream only
    deltaO2Dynamic      =   signal.ts(i)-delayLine(indexO2)*weightsO2';     % dynamic O2 delay in side stream only
    
    delayMMssDynamic    =   deltaMMssDynamic+delayCO2Static;               	% dynamic MMss delay to flow in s
    delayO2Dynamic      =   deltaO2Dynamic+delayCO2Static;                  % dynamic O2 delay to flow in s

    dFlow               =	0;                                              % no delay, flow is reference
    dCO2                =   round(delayCO2Static/dt);                       % static delay CO2 to flow in steps
    dO2                 =   round(delayO2Dynamic/dt);                       % dynamic delay O2 to flow in steps
    dMMss               =   round(delayMMssDynamic/dt);                     % dynamic delay MMss to flow in steps
    
    signal.IvDelay(i) 	=   signal.Iv(i);                                   % copy of signal as flow is reference 
    signal.CO2Delay(i)  =   signal.CO2(min(nMax,i+dCO2));                   % advancing CO2 by dCO2
    signal.O2Delay(i)   =   signal.O2(min(nMax,i+dO2));                     % advancing O2 by dO2
    signal.MMssDelay(i) =   signal.MMss(min(nMax,i+dMMss));                 % advancing MMss by dMMss

    signal.delayMMss(i) =   delayMMssDynamic;                             	% dynamic MMss delay to flow signal
    signal.delayO2(i)   =   delayO2Dynamic;                               	% dynamic O2 delay to flow signal
end