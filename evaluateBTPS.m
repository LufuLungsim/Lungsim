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
function [parameters] = evaluateBTPS(parameters)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup of constants
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    operation       =   parameters.Operation;
    
    btps            =   operation.BTPS;                     % btps correction activated
    
    if btps
        Tkelvin         =   parameters.Physical.tKelvin;
        
        Tinsp           =   operation.TinsGas;
        Tbody           =   operation.Tbody;
        TexpEff         =   operation.TexpEff;
        
        Pref            =   operation.Pref;
        
        xH2OInsp        =   vaporPressure(Tinsp  )/Pref*operation.hRelInsGas;
        xH2OExp         =   vaporPressure(Tbody  )/Pref*operation.hRelBody;
        xH2OSensor      =   vaporPressure(TexpEff)/Pref;
        
        factorTempInsp  =   (Tbody+Tkelvin)/(Tinsp+  Tkelvin);
        factorTempExp   =   (Tbody+Tkelvin)/(TexpEff+Tkelvin);
        factorPinsp     =   (1-xH2OInsp  )/(1-xH2OExp);
        factorPexp      =   (1-xH2OSensor)/(1-xH2OExp);
        
        factorBTPSinsp  =   factorTempInsp*factorPinsp;
        factorBTPSexp   =   factorTempExp *factorPexp;
    else % no BTPS correction
        factorBTPSinsp  =   1;
        factorBTPSexp   =   1;
    end
    
    parameters.Calibration.factorBTPSinsp   =   factorBTPSinsp;
    parameters.Calibration.factorBTPSexp    =   factorBTPSexp;
end
