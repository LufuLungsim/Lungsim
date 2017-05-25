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
function [signal]=signalCalc(signal,parameters)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % setup of general constants
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    operation           =   parameters.Operation;
    
    btps                =   operation.BTPS;                     % btps correction activated
    dynamic             =   parameters.Simulation.dynamicState;	% static correction
    verb                =   parameters.Simulation.verb;         % verbosity state
    
    N2toArN2            =   operation.xAir(5)/(operation.xAir(5)+operation.xAir(7));
    ArtoArN2            =   1.0-N2toArN2;
    
    Tkelvin             =   parameters.Physical.tKelvin;
    
    factorBTPSinsp      =   parameters.Calibration.factorBTPSinsp;
    factorBTPSexp       =   parameters.Calibration.factorBTPSexp;
    
    torontoFile         =   parameters.Simulation.torontoFile;
    
    aFile               =   parameters.Simulation.aFile;
    bFile               =   parameters.Simulation.bFile;
    dtgmbw              =   parameters.Simulation.dtgmbwAnalysis;
    spirowareType       =   parameters.Operation.spirowareType;
    correctionCO2Dry    =   parameters.Operation.CO22CO2RAW;
    gainCO2Dyn          =   parameters.Operation.gainCO2Dyn;
    offsetCO2Dyn        =   parameters.Operation.offsetCO2Dyn; %%%% *1.055;
    MissPlexy           =   parameters.Simulation.MissPlexy;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % data for dynamic control
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tbody               =   operation.Tbody;
    Tnafion             =   operation.Tnafion;
    Tsensor             =   operation.Tsensor;
    
    Pref                =   operation.Pref;
    
    H2OBody             =   vaporPressure(Tbody  )/Pref*operation.hRelBody;
    H2ONafion           =   vaporPressure(Tnafion)/Pref*operation.hRelEnv;
    H2OSensor           =   vaporPressure(Tsensor)/Pref*operation.hRelBody;
    
    %Determine Gasmixture through algorithmic methods
    if ~(torontoFile||bFile)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Signal calculation for Spiroware type data (nothing necessary for toronto type data)
        % estimation of H2O signal at POI1 (in equlilibrium with local temperature)
        % It is assumed that the inspired air is humidified according to the
        % simulated temperature at POI1 due to the wetted tube walls!
        % estimation of H2O signal at POI3 (in equilibrium with labor conditions due to Nafion tube)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % TODO Check if inspired air is saturated with H2O (it's not)
        % TODO set inspirational H2O concentration to zero weil Wandluft
        signal.H2OPoi1      =   vaporPressure(signal.delayTemp)/Pref;
        signal.H2OPoi3      =	ones(size(signal.ts))*H2ONafion;	% side stream H2O molar fraction (temperature and pressure dependent)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Assumptions:
        % O2 is calibrated such that it shows 100% for gas=O2 even if measured with H2ONafion
        % CO2 is calibrated such that is shows dry value even if measured at BTPS at POI1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        signal.corrPoi3     =   (1-signal.H2OPoi3);                 % reduced fraction due to H2O presence!
        if aFile    % spiroware data type
            signal.corrCO2  =   ones(size(signal.CO2Delay));
            if spirowareType == 2	% corrected CO2 signal (solution before 29.03.12): with this factors, LungSim-Bfile CO2 corresponds to Spiroware-Bfile
                signal.corrCO2  =   signal.corrCO2./((1-H2OSensor)*correctionCO2Dry);	% reset dry correction, i.e., set raw CO2 signal
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % The correction by the factor 1./(1-signal.H2OPoi1) corrects to dry condition
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %        signal.corrCO2	=   signal.corrCO2.*(gainCO2Dyn*signal.O2Filter+offsetCO2Dyn);  % dynamic filtering (Spiroware-style)
            signal.corrCO2	=   signal.corrCO2.*(gainCO2Dyn*signal.O2Delay+offsetCO2Dyn)./(1-signal.H2OPoi1);  % dynamic filtering (Spiroware-style) with water content correction
        else        % former solution for non-spiroware CO2 signals
            signal.corrCO2      =   (1-H2OSensor)./(1-signal.H2OPoi1);	% corrects for "true" H2O fraction at POI1
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % evaluation of signals at POI3 (correction for H2O content)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        signal.MMssPoi3     =   signal.MMssDelay;                  % no correction necessary
        signal.O2Poi3       =   signal.O2Delay.*signal.corrPoi3;	% corrects for H2O fraction at POI3
        signal.CO2Poi3      =	signal.CO2Delay.*signal.corrPoi3.*signal.corrCO2;	% corrects for H2O fraction at POI3 and modification of CO2 sensor calibtation issues
        InertGasPoi3        =   1-signal.O2Poi3-signal.CO2Poi3-signal.H2OPoi3;   % deduced inert molar fraction
        %
        % Gas_Species          =	{'O_2', 'CO_2', 'He', 'SF_6', 'N_2', 'H2O', 'Ar'}
        %
        
        if dtgmbw 
           [signal.N2Poi3,  ...
            signal.ArPoi3,  ...
            signal.SF6]     =   getInertGasConcentration(signal, InertGasPoi3, parameters);

            speciesIndex    =   [1,2,4,5,6,7];
            species         =   [signal.O2Poi3,signal.CO2Poi3,signal.SF6,signal.N2Poi3,signal.H2OPoi3,signal.ArPoi3];
            signal          =   getCalculatedSignal(signal, speciesIndex, species, parameters);
            signal.SF6      =   signal.SF6./signal.corrPoi3;
            signal.SF6Poi3  =   signal.SF6;
        else     
            
            signal.N2Poi3   =   InertGasPoi3*N2toArN2;             % distribution to N2 part
            signal.ArPoi3   =   InertGasPoi3*ArtoArN2;             % distribution to Ar part
        
            speciesIndex    =   [1,2,5,6,7];
            species         =   [signal.O2Poi3,signal.CO2Poi3,signal.N2Poi3,signal.H2OPoi3,signal.ArPoi3];
            signal          =   getCalculatedSignal(signal, speciesIndex, species, parameters);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % BTPS correction
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~bFile
        signal.IvEffBTPS    =   signal.IvDelay+signal.Ivss;    % net flow at POI1 (from body point of view)
    end
    
    %%TODO this is a workaround 
    if bFile      
        signal.ArPoi1=zeros(size(signal.CO2));
        signal.ArPoi3=zeros(size(signal.CO2));
        signal.CO2Poi1=signal.CO2;
        signal.CO2Poi3=signal.CO2;
        %signal.H2OPoi1=signal.H2O;
        %signal.H2OPoi3=signal.H2O;
        signal.N2Poi1=signal.N2;
        signal.N2Poi3=signal.N2;
        signal.O2Poi1=signal.O2;
        signal.O2Poi3=signal.O2;
        if dtgmbw
            signal.SF6Poi1=signal.SF6;
            signal.SF6Poi3=signal.SF6;   
        end
        signal.MMssPoi1=signal.MMss; 
        signal.MMssPoi3=signal.MMss;
        signal.MMssCalc=signal.MMss;
        signal.DiffMMss=signal.MMss;
        signal.IvEffBTPS=signal.Iv;
    end
    
    if btps && ~(bFile || MissPlexy)
        if dynamic
            factorT         =	(Tbody+Tkelvin)./(signal.delayTemp+Tkelvin);
            factorP         =   (1-signal.H2OPoi1)/(1-H2OBody);
            
            signal.IvEffBTPS=   signal.IvEffBTPS.*factorT.*factorP;
        else
            signal.IvEffBTPS=   signal.IvEffBTPS.*((signal.IvEffBTPS>0).*factorBTPSinsp+(signal.IvEffBTPS<=0).*factorBTPSexp);
        end
    end
    
    if btps && verb
        if dynamic
            factorStatic1   =   signal.delayLineFilter*factorBTPSinsp+(1-signal.delayLineFilter)*factorBTPSexp;
            factorStatic2   =   (signal.IvFilter>0)*factorBTPSinsp+(signal.IvFilter<=0)*factorBTPSexp;
            
            figure(601);
            subplot(3,1,1);
            plot(signal.ts,factorT.*factorP,signal.ts,factorStatic1,'k-.',signal.ts,factorStatic2,'k--')
            title('Derived dynamic BTPS factors (delayed and static value for comparison)')
            xlabel('t / s')
            ylabel('factors / -')
            legend('BTPS-factor','BTPS delayed','BTPS static')
            axis([-inf,inf,1.0,1.2]);
            
            subplot(3,1,2);
            plot(signal.ts,[factorT,factorP])
            title('Temperature and saturation based BTPS factors')
            xlabel('t / s')
            ylabel('factors / -')
            legend('T-factor','P-factor')
            axis([-inf,inf,1.0,1.2]);
            
            subplot(3,1,3);
            plot(signal.ts,[signal.delayTempFilter,10*signal.delayLineFilter])
            title('Delayed temperature change and delay line')
            xlabel('t / s')
            ylabel('T / °C, delay / -')
            legend('delay Temp','delay line')
            axis([-inf,inf,0,40]);
        else
            factorStatic1   =   signal.delayLineFilter*factorBTPSinsp+(1-signal.delayLineFilter)*factorBTPSexp;
            factorStatic2   =   (signal.IvFilter>0)*factorBTPSinsp+(signal.IvFilter<=0)*factorBTPSexp;
            
            figure(601);
            subplot(2,1,1);
            plot(signal.ts,factorStatic1,'k-.',signal.ts,factorStatic2,'k--')
            title('delayed and static BTPS factors')
            xlabel('t / s')
            ylabel('factors / -')
            legend('BTPS delayed','BTPS static')
            
            subplot(2,1,2);
            plot(signal.ts,10*signal.delayLineFilter)
            title('Delay line')
            xlabel('t / s')
            ylabel('delay / -')
            legend('delay line')
            axis([-inf,inf,0,15]);
        end
    end
end
    
function signal = getCalculatedSignal(signal, speciesIndex, species, parameters)
        
        kappaAir            =   parameters.Physical.kappaAir;
        Cp                  =   parameters.Physical.Cp;
        Cv                  =   parameters.Physical.Cv;
        Mmol                =   parameters.Physical.Mmol; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % evaluation of MMssCalc and MMssDiff signal at POI3
        %
        % implementation according to Speed of sound measurements in
        % gas-mixtures at varying composition using an ultrasonic gas flow
        % meter with silicon based transducers
        % T. Löfquvist, K. Sokas, J. Delsing, EISLAB, Department of computer
        % science and electrical engineering, Luela Unisversity of Technology,
        % Sweden
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        kappaRatio          =   ((species*Cp(speciesIndex))./(species*Cv(speciesIndex)))/(kappaAir*1);
        signal.MMssCalc     =   (species*Mmol(speciesIndex))./kappaRatio;
        
        signal.DiffMMss     =   signal.MMssPoi3-signal.MMssCalc;    % difference (mostly relevant for DTG)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % evaluation of signals at POI1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        signal.O2Poi1       =   signal.O2Delay;                    %
        signal.CO2Poi1      =	signal.CO2Delay.*signal.corrCO2;	% modifies CO2 sensor calibtation issues
        signal.N2Poi1       =   signal.N2Poi3./signal.corrPoi3;     % distribution to N2 part
        signal.ArPoi1       =   signal.ArPoi3./signal.corrPoi3;     % distribution to Ar part
end

function [N2Poi3, ArPoi3, SF6] = getInertGasConcentration(signal, InertGasPoi3, parameters)
    %Verhältniss N2/Ar ist konstant --> kann als ein gas angesehen werden
    %wodurch die mol masse des inerten gases nur noch durch 2 gase bestimmt
    %wird was sich lösen lässt-->(1) mmss = molarMass SF6 + molarMass (N2+Ar)  
    % (2) Anteil Inerte Gase = Anteil SF6 + Anteil (N2 + Ar)
    
    N2toAr              =   parameters.Operation.xAir(5)/(parameters.Operation.xAir(5)+parameters.Operation.xAir(7));
    ArtoN2              =   1.0-N2toAr;   
    kappaAir            =   parameters.Physical.kappaAir;
    Cp                  =   parameters.Physical.Cp;
    Cv                  =   parameters.Physical.Cv;
    Mmol                =   parameters.Physical.Mmol; 
    
    speciesIndex        =   [1,2,6];
    species             =   [signal.O2Poi3,signal.CO2Poi3,signal.H2OPoi3];
    kappaRatio          =   ((species*Cp(speciesIndex))./(species*Cv(speciesIndex)))/(kappaAir*1);
    signal.MMssCalc     =   (species*Mmol(speciesIndex))./kappaRatio;
    signal.DiffMMss     =   signal.MMssPoi3-signal.MMssCalc;
    
    % Gas_Species       =	{'O_2', 'CO_2', 'He', 'SF_6', 'N_2', 'H2O', 'Ar'}
    CpRatio(1:3)        =   [Cp(4), Cp(5), Cp(7)];
    CvRatio(1:3)        =   [Cv(4), Cv(5), Cp(7)];
    MMolRatio(1:3)      =   [Mmol(4), Mmol(5), Mmol(7)];
    
    kappaRatio          =   (CpRatio./CvRatio)/(kappaAir*1);
  
    MMssRatio           =   (MMolRatio)./kappaRatio;
    MMssRatio(2)        =    MMssRatio(2)*N2toAr+MMssRatio(3)*ArtoN2;
    MMssRatio(3)        =   [];
  
    % b = signal vector     (Missing Molar Mass Signal once O2, CO2 and H20 
    %     are subtracted from the whole signal, InertGasProportion of the 
    %     signal) 
    % A = Matrix with conditions to meet, (Molar Mass and Concentration)
    % x = Resulting vector containing the concentrations of SF6, N2, Ar
    % TODO find a way to force positive concentration values. verify SF6 calculation, it seems noisy
    b = [signal.DiffMMss'; InertGasPoi3'];
    A = [MMssRatio;  [1 1]];

    x = A\b;        %Fractions of the gases in the inert gas. TODO, find a method to force the results to be >= 0
    
    SF6 = x(1,:)';
    N2Poi3 = x(2,:)'*N2toAr;
    ArPoi3 = x(2,:)'*ArtoN2; 
    
        
%% 
%     Calculate the gas concentration while ensuring that the concentration is non negative, does not really work and 
%     takes way to long (ca. 240s)
%
%     SF6 = zeros(size(signal.DiffMMss));
%     N2  = zeros(size(signal.DiffMMss));
%     Ar  = zeros(size(signal.DiffMMss));
%     
%     opts = optimoptions('lsqlin', 'Diagnostics', 'off', 'Display', 'off');
%     tic
%     for i = 1:length(signal.DiffMMss)
%         x = lsqlin(A, [signal.DiffMMss(i) InertGasPoi3(i)], [], [], [], [], [0 0], [0.05 0.85], [], opts);
%         SF6(i) = x(1);
%         N2(i)  = x(2)*N2toAr;
%         Ar(i)  = x(2)*ArtoN2; 
%     end
%     toc
%%
    
end
















