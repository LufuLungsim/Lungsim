%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Artificial advancing of main stream signals to apply filter approach
% also for negative delays (necessary for manually synchronized exports
% from wBreath)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% 3.0, 26. Juli 2012
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [signal]      =   signalAdvance(signal,parameters)

    nSamples    =   parameters.Simulation.nAdvancingSamples;
    
    signal.Iv	=	[signal.Iv(nSamples+1:end);signal.Iv(end)*ones(size(signal.Iv(1:nSamples)))];
    signal.CO2	=	[signal.CO2(nSamples+1:end);signal.CO2(end)*ones(size(signal.CO2(1:nSamples)))];
    signal.MMms	=	[signal.MMms(nSamples+1:end);signal.MMms(end)*ones(size(signal.MMms(1:nSamples)))]; 
    
end