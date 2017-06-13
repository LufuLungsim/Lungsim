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
% 2.4, 20. Oktober 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Parameters]=setParametersSpirowareMassSpec(parameterFile,fromFile)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Operation parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Operation.filterOrder   =    0;             % butterworth filter order (upon data reading)
    Operation.filterCutoff  =   10;             % butterworth filter cutoff frequency
    
    Operation.sampleFlow    =   200.0e-6/60;    % nominal sample flow [m^3/s], used if not in data set
    
    Operation.MMair         =   28.97e-3;       % nominal air MM, used to derive %-based CO2 signal from MM-based value
    Operation.MM2Percent    =   12.01e-3;       % proportionality constant between %-based and MM-based CO2 variations
    
    Operation.o2Air         =    0.21;          % nominal O2 concentration in air (used for He MBW data)
    
    Operation.BTPS          =    1;             % optional BTPS correction upon data reading
    Operation.recalibration =    0;             % optional recalibration correction
    Operation.response      =    1;             % optional response (CO2) correction
    
    Operation.maximalLag    =    0.70;          % maximal expected lag of slowes signal
    Operation.Pref          =	 1.01325e5;     % normal pressure in Pa
    Operation.Tenv          =   25.0;           % environmental temperature in °C
    Operation.hRelEnv       =	 0.50;          % reference relative humidity in lab air (surrounding Nafion tube)
    Operation.TinsGas       =  Operation.Tenv;  % reference temperature of inspiration gas in °C
    Operation.hRelInsGas    =	 0.03;          % reference relative humidity of inspiration gas
    
    Operation.Tbody         =   37.0;           % body temperature
    Operation.hRelBody      =	 1.0;           % reference relative humidity during expiration
    
    Operation.TexpEff       =   33.0;           % eff. temperature of air during expiration at sensor
    Operation.Tsensor       =   28.0;           % sensor reference temperature
    Operation.Tbody1        =   Operation.Tbody;% body reference temperature 1
    Operation.Tbody2        =   33.0;           % body reference temperature 2
    
    Operation.Tnafion       =	Operation.Tenv;	% reference temperature of Nafion tube in °C
    Operation.xAir          =   [0.20943,0.00039,0.0,0.0,0.78084,0,0.00934]'; % composition of air
    Operation.xCal          =   [0.15000,0.05000,0.0,0.0,0.80000,0,0.00000]'; % composition of calibration gas
    
    Operation.toronto       =   [0.00000,0.00000,0.04,0.04,0.00000,0.00000,0.00000]'; % composition of air (Attention: only He and SF6 is set!!)
    
    Operation.dateCO2CORR   =   datenum('29-Mar-2012'); % before this date, data was recoreded as CO2 (corrected for dry conditions)
    Operation.dateCO2RAW    =   datenum('9-May-2012');  % from this date, data was recoreded as CO2_RAW (uncorrected, i.e., raw sensor value)
    Operation.CO22CO2RAW    =   1.075;          % transformation of CO2 signal back to CO2_RAW signal (CO2 under dry conditions)
    Operation.gainCO2Dyn    =   0.0009*100;     % Spiroware value is 0.0009 as molar fractions are in % (not ratios)
    Operation.offsetCO2Dyn  =   0.9855;         % The formula correctionFactor = gainCO2Dyn * O2 + offsetCO2Dyn equals 1 for O2=0.1611 (expiration)
    
    Operation.spirowareType =  -1;              % 0: after Operation.dateCO2RAW, 2: before Operation.dateCO2CORR, 1: in between these dates, -1: no Spiroware type
    
    Operation.nVentilation  =  10;              % #breaths for minute ventilation evaluation
    
    Operation.deltaCO2O2    =   0;              % manual correction of delay CO2->O2
    Operation.deltaCO2MMss  =   0;              % manual correction of delay CO2->MMss
    Operation.dTauCO2       =   0;              % manual correction of tau CO2
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Physical parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Physical.tKelvin        =    273.15;        % absolut zero on Kelvin scale / °C
    Physical.k              =    0.025;         % heat conduction constant / W/K/m
    Physical.D              =    0.2e-4;        % CO2/N2-O2 diffusion constant / m/s**2
    Physical.rho            =    1.18;          % density air / kg/m**3
    Physical.Rgas           =    8.314;         % universal gas constant / J/K/mol
    Physical.Mair           =    28.8e-3;       % kg/mol
    Physical.cpRho          = 	 7/2*Physical.Rgas/Physical.Mair*Physical.rho; % heat capacity / J/K/kg
    Physical.c0             =	 Operation.Pref/(Physical.Rgas*(Operation.Tenv+273.0)); % mean molar concentration
    
    Physical.Names          =	{'O_2', 'CO_2', 'He', 'SF_6', 'N_2', 'H2O', 'Ar'};
    Physical.Mmol           =	[32.0 44.01 4.0 146.0 28.013 18.015 39.948]'*1e-3;      % molar masses of the species kg/mol
    %Physical.kappa          =	[1.3964 1.2976 1.63 1.0935 1.4014 1.33 1.6757]';        % adiabatic constant (value for water from wikipedia.en)
    Physical.kappaAir       =   1.4012;                                                 % adiabatic constant of air (consistent with xAir, Cp and Cv)
    
    Physical.Cp             =   [29.38 36.94 20.79 97.85 29.12 37.47 20.79]';           % ratio, relevant for Mstar calculation (p. 83 Diss Buess)
    Physical.Cv             =   [21.00 28.48 12.47 89.54 20.80 28.03 12.47]';           % ratio, relevant for Mstar calculation (p. 83 Diss Buess)
    Physical.R              =   [259.827 188.922 2077.058 8.314472 296.839 461.401 208.122];
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Simulation related parameters: discretization settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Simulation.fullOutput               =   1;              % output switch
    Simulation.reducedOutput            =   0;              % output switch
    
    Simulation.verb                     =    0;             % verbosity switch (controls degree of debugging output)
    Simulation.compatibility            =    0;             % compatibility switch (sets pre 4.8.1 variant of MBW end detection in MassSpec)
    Simulation.dt                       =    0.005;         % time increment (used only for visualization)
    Simulation.nAdvancingSamples        =    0;             % advancing of Iv and CO2 in case of syncronization overshooting
    
    Simulation.bfileOutput              =    1;             % output B-file style
    Simulation.cfileOutput              =    0;             % output C-file style
    Simulation.inselOutput              =    0;             % output I-file style
    Simulation.leicesterOutput          =    0;             % output L-file style
    
    Simulation.slowState                =    0;             % slowing of signals
    Simulation.fastState                =    1;             % accelerating of signals
    
    Simulation.torontoFile              =    1;             % toronto file (related to mass spectrometer data (default at the moment)
    
    Simulation.aFile                    =    1;             % new A file (with 'A-....' file name)
    Simulation.bFile                    =    0;             % new A file (with 'A-....' file name)
    Simulation.sbwMm                    =    0;             % old SBW file (CO2 in g/mol)
    Simulation.sbwPercent               =    0;             % old SBW file (CO2 in %)
    Simulation.heMbw                    =    0;             % old He-based MBW file (CO2 in %)
    
    Simulation.graphState               =    0;             % no graphs are produced during batch operations
    Simulation.interactiveYes           =    0;             % allow for interactive slope fits
    Simulation.interactiveNo            =    1;             % allow for interactive slope fits
    Simulation.interactiveDeviation     =    0;             % allow for interactive slope fits
    Simulation.rmsCrit                  =    0.0001;        % critical rms value to trigger interactive processing
    Simulation.tauFix                   =    0.04;
    
    Simulation.fileCalibrationNo        =    1;             % per file calibration
    Simulation.fileCalibrationYes       =    0;             % per file calibration
    Simulation.fileCalibrationModifier	=    '-;_pre';      % per file calibration modifier
    
    Simulation.n2mbwAnalysis            =    1;             % n2mbw analysis type
    Simulation.dtgsbwAnalysis           =    0;             % dtgsbw analysis type
    Simulation.dtgmbwAnalysis           =    0;             % dtgmbw analysis type
    Simulation.dtgmbw                   =    0;             % Double Tracer Gas Multiple Breath Washout
    Simulation.dtgmbwBaby               =    0;             % Double Tracer Gas Multiple Breath Washout
    Simulation.fastSlow                 =    0;             % Analyse Fast Slow compartments
    Simulation.fastSlowPredict          =    0;             % Predict Fast Slow compartments
    Simulation.plainAnalysis            =    0;             % plain analyis type
    Simulation.conversionAnalysis       =    0;             % data conversion type
    Simulation.tidalAnalysis            =    0;             % tidal analysis
    Simulation.MissPlexy                =    0;             % If the measurements were done in Miss Plexy
    Simulation.Set1                     =    0;             % If the measurements were done with Set1
    Simulation.Set2                     =    1;             % If the measurements were done with Set2 or 3
    
    Simulation.staticState              =    0;             % static vs. dynamic BTPS correction
    Simulation.dynamicState             =    1;             % static vs. dynamic BTPS correction
    
    Simulation.slopesState              =    0;             % include slopes in log file
    Simulation.minVentState             =    0;             % include minute ventilations in log file
    
    Simulation.gainCO2MM                =    1/(0.04401-0.02897);	% raw gain CO2MM
    
    Simulation.relativeError            =    2.0e-4;        % relative error during iterations
    Simulation.maximalIteration         =	50;             % maximla iteration count
    
    Simulation.versionString            =   'Version 4.11.8';  % version string
    
    Simulation.lciStandard              =    1;             % Evaluation of LCI by standard procedure
    Simulation.lciByFit                 =    0;             % Evaluation of LCI by fitting N2et(Vexp)
    
    Simulation.ScondSacin               =    1;             % Evaluation of Scond/Sacin activated
    Simulation.exclusionLimit           =    0.25;          % limit of deviation of Vtidal from mean value
    Simulation.onlySmall                =    1;             % exclude only too small tidal breaths
    Simulation.ScondSacinCheck          =    0;             % Interactive check of ScondSacin included breaths after completion of all breaths
    Simulation.ScondSacinNorm           =    0;             % Normalization of slope (0=SnIII*Vexp, 1=SnIII)
    Simulation.lBound                   =    1.5;           % lower bound of TO to evaluated Scond/ScondStar
    Simulation.uBound                   =    6.0;           % upper bound to TO to evaluated Scond
    Simulation.uBoundStar               =    3.0;           % upper bound to TO to evaluated ScondStar
    Simulation.lBoundStar               =    0;             % lower bound to TO to evaluated ScondSta
    
    Simulation.CapnoIndices             =    0;             % Evaluation of capno indices according to "Volumetric Capnography in Infants with BPD", J. pediatr 164, p283-8 (2014)
    Simulation.CapnoCheck               =    0;             % Interactive check of capno indices after completion of all breaths
    
    Simulation.tidalMeanState           =    0;             % include tidal means in log file
    Simulation.fromBreathTidal          =    11;            % start breath to evaluate mean
    Simulation.nBreathTidal             =    10;            % # breaths to evaluate mean
    
    Simulation.TOevaluation             =    0;             % TO based index (alternative to LCI)
    
    Simulation.resultDirPath            =   getenv('TEMP'); % default result dir Path
    
    Simulation.nameHash                 =   cell(0,2);      % empty hash-array with correct structure
    
    Simulation.settingsState            =   0;              % calibration data form settings file
    Simulation.manualState              =   0;              % manual calibration (anew for each signal file)
    Simulation.autoState                =   1;              % automatic evaluation of calibration (not yet implemented)
    
    Simulation.croppingState            =   0;              % signals are manually cropped before processing
    Simulation.onlyPlotState            =   0;              % signals are only plotted
    
    Simulation.nBreathMax               =  50;              % maximal number of recorded breaths in log file
    
    Simulation.directoryCalibration     =   '';             % default directory to start search for calibration file
    Simulation.directoryDatafiles       =   '';             % default directory to start search for data files
    
    Simulation.MMeeMethod               =   'MMeeR98';
    Simulation.MMeeTMin                 =   0.65;
    Simulation.MMeeTMax                 =   0.95;
    
    Simulation.LCIbyFitAutoVal          =   0;
    Simulation.LCIbyFitManVal           =   0;
    Simulation.LCIbyFitMaxStdDev        =   0.05;
    
    Simulation.MomentRatio              =   0;
    Simulation.LinearByFit              =   1;
    Simulation.ExponentialByFit         =   0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Device geometrical and physical parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Device.length       =   [0.06,0.14,1.30];                   % device length per branch / m
    Device.radius       =   [0.0075,0.0075,0.0008];             % tube radii per branche / m
    Device.area         =   (Device.radius.^2)*pi;              % tube section area per branche / m**2
    
    Device.vChamber     =   1.2e-6;                             % side chamber volume / m**3
    
    Device.dFlow        =   2*Device.radius(2)*sqrt(2);         % length of sensing path / m
    Device.dCam         =   [0.02,0.02];                        % length of side chambers sensing path / m
    Device.alpha        =   5.0*1;                              % heat transfer coefficient / W/K/m**2
    Device.gamma        =   1.0e-5;                             % side chamber coupling constant / m**3/s
    
    Device.volumePostcap=   26.9e-6;                            % postcap volume
    Device.volumePrecap =   37.9e-6;                            % precap volume
    Device.volumeBody2  =   45.0e-6;                            % body volume 2
    Device.volumeBody1  =   45.0e-6;                            % body volume 1
    
    Device.tauTemp      =    2.0;                               % precap tau relaxation
    Device.tauBody2     =    0.5;                           	% body tau relaxation 2
    Device.tauBody1     =    0.5;                               % body tau relaxation 1
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calibration data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Calibration.settingsFilePath    =   'not yet set';
    Calibration.filePath            =   'not yet set';
    Calibration.delayFlow           =	0.01;                                  	% delay of flow USFM to species signals (not detected automatically)
    Calibration.delayCO2    =	0.01;           % delay of CO2 to flow (main stream)
    Calibration.delayO2    	=   0.683;          % delay of O2 to flow (main stream)
    Calibration.delayMMss  	=   0.794;          % delay of MMss to flow (main stream)
    Calibration.coefficientsOptimal =   [];
    Calibration.time.min            =   0;
    Calibration.time.max            =   0;
    
    Calibration.endInsO2            =   Operation.xAir(1);                      % end inspiration O2 nominal value
    Calibration.endInsCO2           =   Operation.xAir(2);                      % end inspiration CO2 nominal value
    Calibration.endInsMMss          =   (Operation.xAir'*Physical.Mmol);        % no kappa-correction necessary
    
    Calibration.factorBTPSinsp      =   1.0;                                    % set during operation
    Calibration.factorBTPSexp       =   1.0;                                    % set during operation
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Packing the setup variables into the structure Parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Parameters.Physical     =   Physical;
    Parameters.Device       =   Device;
    Parameters.Operation	=   Operation;
    Parameters.Simulation	=   Simulation;
    Parameters.Calibration  =   Calibration;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Basic settings
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if fromFile
        if exist(parameterFile,'file')
            ParametersLoaded = load(parameterFile,'parameters');
            if ~isempty(Parameters)
                %Check if the loaded Parameters have the same number of
                %fields as the standard values.
                if length(fieldnames(ParametersLoaded)) == length(fieldnames(Parameters))
                    Parameters = ParametersLoaded;
                end
            end
        else
            fprintf('cannot load file %s; using standard values\n',parameterFile);
        end
    end
    
    
end
