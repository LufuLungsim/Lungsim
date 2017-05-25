%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Matlab function to set signal correction parameters
%
% It calculates especially the step response for different convective velocity values
% The step consists simultanuously of a sudden concentration change as well
% as a temperatur step during expiration
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% 2.5, 30. Januar 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Parameters]=setParametersWbreath(parameterFile,fromFile)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Operation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Operation.filterOrder   =    0;             % butterworth filter order (upon data reading)
    Operation.filterCutoff  =   10;             % butterworth filter cutoff frequency
    
    Operation.o2Air         =    0.21;          % nominal O2 concentration in air (used for He MBW data)
    
    Operation.Pref          =	 1.01325e5;     % normal pressure in Pa
    Operation.Tenv          =   25.0;           % environmental temperature in °C
    Operation.TinsGas       =  Operation.Tenv;  % reference temperature of inspiration gas in °C
    Operation.hRelInsGas    =	 0.03;          % reference relative humidity of inspiration gas
    
    Operation.Tbody         =   37.0;           % body temperature
    Operation.hRelBody      =	 1.0;           % reference relative humidity during expiration
    
    Operation.xAir          =   [0.20943,0.00039,0.0,0.0, 0.78084,0,0.00934]'; % composition of air
    Operation.xTG           =   [0.20943,0.00039,0.0,0.04,0.74084,0,0.00934]'; % composition of air
    
    Operation.nVentilation  =  10;              % #breaths for minute ventilation evaluation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Physical parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Physical.tKelvin        =  273.15;          % absolut zero on Kelvin scale / °C
    Physical.k              =    0.025;         % heat conduction constant / W/K/m
    Physical.D              =    0.2e-4;        % CO2/N2-O2 diffusion constant / m/s**2
    Physical.rho            =    1.18;          % density air / kg/m**3
    Physical.Rgas           =    8.314;         % universal gas constant / J/K/mol
    Physical.Mair           =    28.8e-3;       % kg/mol
    Physical.cpRho          = 	7/2*Physical.Rgas/Physical.Mair*Physical.rho; % heat capacity / J/K/kg
    Physical.c0             =	Operation.Pref/(Physical.Rgas*(Operation.Tenv+273.0)); % mean molar concentration
    
    Physical.Names          =	{'O_2', 'CO_2', 'He', 'SF_6', 'N_2', 'H2O', 'Ar'};
    Physical.Mmol           =	[32.0,44.01,4.0,146.0,28.013,18.015,39.948]'*1e-3;      % molar masses of the species kg/mol
    %Physical.kappa          =	[1.3964 1.2976 1.63 1.0935 1.4014 1.33 1.6757]';        % adiabatic constant (value for water from wikipedia.en)
    Physical.kappaAir       =   1.4012;                                                 % adiabatic constant of air (consistent with xAir, Cp and Cv)
    
    Physical.Cp             =   [29.38,36.94,20.79,97.85,29.12,37.47,20.79]';           % ratio, relevant for Mstar calculation (p. 83 Diss Buess)
    Physical.Cv             =   [21.00,28.48,12.47,89.54,20.80,28.03,12.47]';           % ratio, relevant for Mstar calculation (p. 83 Diss Buess)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Operation parameters (additional)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Operation.MMair         =   Physical.Mmol'*Operation.xAir;                          % nominal air MM
    Operation.MMTG          =   Physical.Mmol'*Operation.xTG;                           % nominal tracer gas MM
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Device.volumePrecap                 =    4.0e-6;     	% precap volume
    Device.volumePostcap                =    0.9e-6;      	% postcap volume
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation related parameters: discretization settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Simulation.verb                     =    0;             % verbosity switch (controls degree of debugging output
    Simulation.dt                       =    0.005;         % time increment (used only for visualization)
    
    Simulation.graphState               =    1;             % no graphs are produced during batch operations
    Simulation.interactiveYes           =    0;             % allow for interactive slope fits
    Simulation.interactiveNo            =    1;             % allow for interactive slope fits
    Simulation.interactiveDeviation     =    0;             % allow for interactive slope fits
    Simulation.rmsCrit                  =    0.0001;        % critical rms value to trigger interactive processing
    Simulation.tauFix                   =    0.04;
    
    Simulation.slopesState              =    0;             % include slopes in log file
    Simulation.minVentState             =    0;             % include minute ventilations in log file
    
    Simulation.versionString            =   'Version 2.1.2';  % version string
    
    Simulation.SF6                      =    1;             % MM based on SF6=1, He=0
    Simulation.sideStream               =    0;             % side stream data used
    
    Simulation.lciStandard              =    1;             % Evaluation of LCI by standard procedure
    Simulation.lciByFit                 =    0;             % Evaluation of LCI by fitting N2et(Vexp)
    Simulation.TOevaluation             =    0;             % TO based index (alternative to LCI)
    Simulation.momentRatios             =    1;             % Moment ratio evaluation
    
    Simulation.washout                  =    1;             % washout state
    
    Simulation.meanEE                   =    1;             % mean=1 or median=2 or byFit=3 to evaluate ee values
    Simulation.meanEI                   =    1;             % mean=1 or median=2 or byFit=3 to evaluate ei values
    Simulation.lowerBoundEE             =    0.65;          % lower bound (volumina based)
    Simulation.upperBoundEE             =    0.95;          % upper bound (volumina based)
    Simulation.lowerBoundEI             =    0.45;          % lower bound (volumina based)
    Simulation.upperBoundEI             =    0.75;          % upper bound (volumina based)
    
    Simulation.deltaMMState             =    0;             % deltaMMState=1: using given deltaMMFixed value instead of data derived
    Simulation.deltaMMFixed             =    0.0056;        % default of fixed deltaMM value
    
    Simulation.tidalMeanState           =    0;             % include tidal means in log file
    Simulation.fromBreathTidal          =    11;            % start breath to evaluate mean
    Simulation.nBreathTidal             =    10;            % # breaths to evaluate mean
    
    Simulation.ScondSacin               =    1;             % Evaluation of Scond/Sacin activated
    Simulation.exclusionLimit           =    0.25;          % limit of deviation of Vtidal from mean value
    Simulation.onlySmall                =    1;             % exclude only too small tidal breaths
    Simulation.ScondSacinCheck          =    0;             % Interactive check of ScondSacin included breaths after completion of all breaths
    Simulation.lBound                   =    1.5;           % lower bound of TO to evaluated Scond/ScondStar
    Simulation.uBound                   =    6.0;           % upper bound to TO to evaluated Scond
    Simulation.uBoundStar               =    3.0;           % upper bound to TO to evaluated ScondStar
    Simulation.lBoundStar               =    0;             % lower bound to TO to evaluated ScondStar
    
    Simulation.resultDirPath            =   getenv('TEMP'); % default result dir Path
    
    Simulation.nameHash                 =   cell(0,2);      % empty hash-array with correct structure
    
    Simulation.settingsState            =   1;              % calibration data form settings file
    Simulation.manualState              =   0;              % manual calibration (anew for each signal file)
    Simulation.autoState                =   0;              % automatic evaluation of calibration (not yet implemented)
    
    Simulation.croppingState            =   0;              % signals are manually cropped before processing
    Simulation.tidalState               =   0;              % signals are treated as tidal maneuver
    Simulation.onlyPlotState            =   0;              % signals are only plotted
    
    Simulation.nBreathMax               =  50;              % maximal number of recorded breaths in log file
    
    Simulation.directoryCalibration     =   '';             % default directory to start search for calibration file
    Simulation.directoryDatafiles       =   '';             % default directory to start search for data files
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Calibration.settingsFilePath    =   'not yet set';
    Calibration.filePath            =   'not yet set';
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Packing the setup variables into the structure Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Parameters.Physical     =   Physical;
    Parameters.Operation	=   Operation;
    Parameters.Device   	=   Device;
    Parameters.Simulation	=   Simulation;
    Parameters.Calibration  =   Calibration;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Basic settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if fromFile
        if exist(parameterFile,'file')
            ParametersLoaded = load(parameterFile,'parameters');
            if ~isempty(Parameters)
                if length(fieldnames(ParametersLoaded)) == length(fieldnames(Parameters))
                    Parameters = ParametersLoaded;
                end
            end
        else
            fprintf('cannot load file %s; using standard values\n',parameterFile);
        end
    end
    
end


