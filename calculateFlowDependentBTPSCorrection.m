%

function [  ] = calculateFlowDependentBTPSCorrection( signal,parameters)
   
    dtgmbw              = parameters.Simulation.dtgmbwAnalysis;
    %Set Physical Parameters
    % Gas_Species       = {'O_2', 'CO_2', 'He', 'SF_6', 'N_2', 'H2O', 'Ar'}
    if dtgmbw == 1
        speciesReact	= [1,2,6];
        speciesInert    = [4,5,7];
    else
        speciesReact	= [1,2];
        speciesInert    = [6,5,7];
    end
    species             = [speciesReact, speciesInert];
    T0                  = parameters.Physical.tKelvin;      % 0°C in Kelvin
    Tamb                = 20;                               % Ambient Temperature in °C
    Tsatp               = T0+Tamb;
    p0                  = parameters.Operation.Pref;      
    
    %Correct the measured signal
    [signal]            = dataDelay(signal,parameters);  
        
    %Set and initialize Simulation Parameters
    refMMss             = true;
    iterations          = 5;
    signal.H2O          = ones(size(signal.O2Filter)).*vaporPressure(Tamb)/p0;     %Initial guess
    err                 = zeros(iterations*2, length(signal.O2Filter));
    TCorrection         = zeros(size(signal.MMms));
   
    for i = 1:iterations
        
        %Correct for Gas content (signal.H2O)
        [errMM, errC, signal] = calculateMolarMass(signal, TCorrection, speciesReact, speciesInert, Tsatp, parameters, refMMss);

        err(i*2-1,:) = errMM;
        %signal.H2O = min(max(signal.H2O+errC,0), vaporPressure(Tamb-TCorrection)/p0);  
       
        %Correct for Temperature
        [errMM, errC, signal] = calculateMolarMass(signal, TCorrection, speciesReact, speciesInert, Tamb, parameters, refMMss);
        
        cSpecies    = getSpeciesMatrix(signal, dtgmbw);
        MMcalc      = calculateMM(species, cSpecies, TCorrection, Tsatp, parameters);

        if refMMss
            err(i*2,:)  = errMM;
            Terr        = Tsatp*(MMcalc.^2./(MMcalc.*signal.MMssFilter) - 1); 
        else
            err(i*2,:)  = errMM;
            Terr        = Tsatp*(MMcalc.^2./(MMcalc.*signal.MMms) - 1);   
        end   
        TCorrection = TCorrection+Terr;
        %TCorrection = zeros(size(signal.MMms));
       
        plotStep(Terr, TCorrection, cSpecies, signal, dtgmbw, MMcalc, err)
        
        pause(0.05)
    end 
end


% Difference in the calculated Molar Mass and the measured molar Mass 
% as calculated by following formula
% 
% errMM     = MMss - MMcalc * kAir/kCalc
%
% Derived from te formulas:
%
% MMss/kAir = MMcalc/kCalc
% MMcalc    = cSpecies*MMSpecies
% kCalc     = ((cSpecies*kSpecies)/(kSpecies-1)) / (cSpecies/kSpecies-1))
% kCalc is calculated following the formula given by https://de.wikipedia.org/wiki/Gasgemisch
function [errMM, errC, signal] = calculateMolarMass(signal, TCorrection, speciesReact, speciesInert, Tsatp, parameters, refMMss)
                      
    dtgmbw      = parameters.Simulation.dtgmbwAnalysis;
    
    if refMMss
        refMM  = signal.MMssFilter;
    else
        refMM  =  signal.MMms;
    end
    
    if dtgmbw
        gasReact        = [signal.O2Filter, signal.CO2Filter, signal.H2O];
    else
        gasReact        = [signal.O2Filter, signal.CO2Filter];
    end  
    % Given Parameters
    kAir        = parameters.Physical.kappaAir;
    N2ToAr      = parameters.Operation.xAir(5)/(parameters.Operation.xAir(5)+parameters.Operation.xAir(7));
    ArToN2      = 1-N2ToAr;
    
    % Gas_Species       =	{'O_2', 'CO_2', 'He', 'SF_6', 'N_2', 'H2O', 'Ar'}
    % Species parameters
    inertGasRatio = ones(size(signal.O2Filter))-sum(gasReact, 2);
    Cp          = parameters.Physical.Cp';
    Cv          = parameters.Physical.Cv';
    kSpecies    = Cp([speciesReact, speciesInert])./Cv([speciesReact, speciesInert]);                       % Kappa Species
    MMSpecies   = parameters.Physical.Mmol';     % Molar Mass Species
    

    CpReact         =  Cp(speciesReact);        
    CvReact         =  Cv(speciesReact);   
    kReact          =  getKappaGasMixture(CpReact./CvReact, gasReact);
    MMReact         =  MMSpecies(speciesReact);
    MMReact         =  sum(bsxfun(@times,gasReact,MMReact), 2);%    bsxfun(@times,gasReact, MMReact);
    inertGasMass    =  refMM-MMReact./kReact.*kAir;

    MMRatio(1:3)    =  MMSpecies(speciesInert);
    
    MMRatio(2)      =  MMRatio(2)*N2ToAr+MMRatio(3)*ArToN2;
    MMRatio(3)      = [ ];

    A               = [MMRatio;  [1 1]];
    b               = [inertGasMass'; inertGasRatio'];       

    x = A\b; 
    
    if dtgmbw    
        signal.SF6      = max(x(1,:),0)';
    else
        signal.H2O      = max(x(1,:),0)';
    end
    signal.N2       = min(max(inertGasRatio.*N2ToAr, 0), 0.81);
    signal.Ar       = min(max(inertGasRatio.*ArToN2, 0), 0.1); 
    
    cSpecies        = getSpeciesMatrix(signal, dtgmbw);
    MMcalc          = calculateMM([speciesReact, speciesInert], cSpecies, TCorrection, Tsatp, parameters);
    kappaCalc       = getKappaGasMixture(kSpecies, cSpecies);
    
    % Error calculation Formula
    errMM   = refMM - MMcalc.*kAir./kappaCalc;    
    errC = 1-sum(cSpecies,2);
end

function kappaCalc = getKappaGasMixture(kSpecies, cSpecies)
    kappaCalc       = sum(bsxfun(@times,1./(kSpecies-1),bsxfun(@times,kSpecies,cSpecies)), 2)./sum(bsxfun(@times,1./(kSpecies-1),cSpecies),2);
end

function MMcalc = calculateMM(species, cSpecies, TCorrection, Tsatp, parameters)
    MMSpecies           = parameters.Physical.Mmol(species)';
   
    % Molar Mass Species
    MMcalc    = sum(bsxfun(@times,cSpecies, MMSpecies), 2);         %signal.H2O correction
    MMcalc    = MMcalc./(1+TCorrection./Tsatp);                     %Temperature Correction
end

function cSpecies = getSpeciesMatrix(signal,dtgmbw)
    if dtgmbw   
        cSpecies    = [ signal.O2Filter,     ...
                        signal.CO2Filter,    ...
                        signal.SF6,         ...
                        signal.N2,          ...
                        signal.H2O,               ... 
                        signal.Ar];
    else        
        cSpecies    = [ signal.O2Filter,     ...
                        signal.CO2Filter,    ...
                        signal.N2,          ...
                        signal.H2O,               ... 
                        signal.Ar];
    end
end

function plotStep(Terr, TCorrection, cSpecies, signal, dtgmbw, MMcalc, err)
        hTemp = getOrMakeFigure('Signal Development');
        
        subplot(2,2,1)  
        hold off
        plot(Terr)
        hold on
        plot(TCorrection)
        legend('Terr', 'Temperature Correction')
        title('Temperature error (Based on Molar Mass difference)')
        ylabel('[°C] needed to correct from MM Measured to MMCalc')
        xlabel('Measurement points (n)')
        
        
        subplot(2,2,2) 
        hold off
        plot(sum(cSpecies,2))
        hold on
        plot(signal.Ar)
        plot(signal.N2)
        plot(signal.O2Filter)
        plot(signal.CO2Filter)
        plot(signal.H2O);
        if dtgmbw
            plot(signal.SF6)
            legend('Sum', 'Ar', 'N2', 'O2', 'CO2', 'H2O', 'SF6')
        else
            legend('Sum', 'Ar', 'N2', 'O2', 'CO2', 'H2O')
        end
        title('Concentration Estimate')
        ylabel('Concentration Ratio [n/1]')
        xlabel('Measurement points (n)')
        
        subplot(2,2,3);
        hold off
        plot(signal.MMssFilter)
        hold on
        plot(signal.MMms)
        plot(MMcalc) 
        legend('MMss', 'MMms', 'MMcalc');
        title('Molar Mass Comparison')
        ylabel('[kg/mol]')
        xlabel('Measurement points (n)')
        
        subplot(2,2,4);
        hold off
        plot(mean(err,2))
        hold on
        plot(std(err,0,2))
        
        legend('Mean Error', 'std Error')  
        title('Error Development \newline \fontsize{10} errMM = MMMeasured - MMcalc.*kAir./kappaCalc')
        ylabel('Error [n]')
        xlabel('Optimisation Iterations (n)')
end
