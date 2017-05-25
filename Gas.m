classdef Gas
    %GAS Summary of this class goes here
    %   CO2             = Gas                   %Initialisierten des Gases
    %   CO2.name        = 'CO2'                 %Zuweisen eines Namens
    %   name            = CO2.name              %Auslesen des Namens
    %   name            = CO2.getName           %Aufruf der getName funktion
    %   name            = getName(CO2)          %Aufruf der getName funktion
    %   CO2.Standard.lci= 1                     %Zuweisen eines wertes in einer Struktur
    %   lciN2           = CO2.Standard.lci      %Auslesen eines wertes in einer Struktur
    
    properties
        name        = 0;                        %Name of the gas
        et          = 0;                        %end tidal value (expirations)
        eti         = 0;                        %end tidal value (inspirations)
        exp         = 0;                        %expired volume (of tracer gas)
        insp        = 0;                        %inspired volume
        expRoof     = 0;                        %expired volume of tracer gas in a washin calculated as normal but flipped
        inspRoof    = 0;                        %inspired volume of tracer gas in a washin calculated as normal but flipped
        meanInsp    = 0;                        %mean inspired volume (used for mass spectrometer based analyses only)
        reInsp      = 0;                        %re-inspired volume
        notReExp    = 0;                        %re-expired volume
        slopePlain  = 0;                        %plain slope during expiration
        slope       = 0;                        %normalized slope during expiration (by et species concentration)
        calcSlope   = 0;
        diffSlope   = 0;
        cetStart    = 0;                        %start value of cet
        normFactor  = 0;                        %Value to which the tracer gas is normed
        capno       = 0;
        breaths     = {};
        fastSlow    = 0;                        %If the fastslow decomposition is made, save the struct to fastSlow
        twoCompartments = cell(1);
        threeCompartments = cell(1);
        General=struct(...                      %diagnostic information without direct relation to the evaluation method (standard or byFit)
            'netExpVol', 0, ...                 %net expired tracer gas volume per breath (was called VolN2Netto)
            'cumNetExpVol', 0, ...              %cummulative expired net tracer gas volume per breath (CumVolN2Netto)
            'cet_start', 0, ...                 %start value of tracer gas for MBW
            'cev', 0, ...                       %cummulative expired volume per breath (all species)
            'cev_ds', 0, ...                    %CEV minus full dead space per breath (pre and post cap)
            'frcsp', 0, ...                     %FRC per breath
            'frcao', 0, ...                     %FRC minus pre cap volume per breath
            'to', 0, ...                        %turnover per breath per breath (relative to classic FRC)
            'logCet_norm', 0, ...              	%logarithm of normed tracer gas end tidal values
            'cet_normFitUpper', 0, ...         	%fit to upper part of end tidal concentration curve
            'cet_normFitLower', 0, ...        	%fit to lower part of end tidal concentration curve
            'minuteVentilation', 0, ...         %ventilation per minute per breath
            'minuteVentilationTotal', 0, ...	%total ventilation per minute per breath
            'minuteVentilationCrit', 0, ...     %critical ventilation per minute per breath
            'breathRateTotal', 0, ...           %total breath rate
            'breathRateCrit', 0, ...            %critical breath rate
            'scondSacin', {});                  %scond/sacin data structure
        
        Standard=struct(...                     %standard evaluation of FRC and LCI for each critical end ratio, including 2.5% value (was called Direct)
            'cet_norm', 0, ...                 	%normed tracer gas end tidal values per critical end ratio
            'frc', 0, ...                      	%FRC per critical end ratio
            'lci', 0, ...                      	%LCI per critical end ratio
            'duration', 0, ...                 	%maneuver duration per critical end ratio
            'minuteVentilation', 0, ...       	%mean ventilation per minute per critical end ratio
            'nBreaths', 0, ...                 	%number of breaths per critical end ratio
            'breathRate', 0, ...               	%mean breath rate per critical end ratio
            'momentRatio', 0);                  %Moment ratio per cutoff

        ByFit=struct(...                        %by Fit evaluation of FRC and LCI for each critical end ratio, including 2.5% value (was called Direct)
            'cet_norm', 0, ...                  %normed tracer gas end tidal values per critical end ratio
            'frc', 0, ...                       %FRC per critical end ratio
            'lci', 0, ...                       %LCI per critical end ratio
            'duration', 0, ...                  %maneuver duration per critical end ratio
            'minuteVentilation', 0, ...         %mean ventilation per minute per critical end ratio
            'partialBreath', 0, ...             %the ratio between the breaths between which the end tidal value was found
            'nBreaths', 0, ...                  %breath index signifying LCI per critical end ratio
            'breathRate', 0, ...                %mean breath rate per critical end ratio
            'deltaVolume', 0, ...               %Difference between CEV standart and CEV by Fit
            'momentRatio', 0);                  %Moment ratio per cutoff
            
        Turnover=struct(...                     %evaluation of et(Gas) for a given TO
            'cet_norm', 0, ...                  %normed tracer gas end tidal values per critical end ratio
            'frc', 0, ...                       %FRC per critical end ratio
            'lci', 0, ...                       %LCI per critical end ratio
            'duration', 0, ...                  %maneuver duration per critical end ratio
            'nBreaths', 0);                     %breath index signifying LCI per critical end ratio
        
    end
    
    methods
        function name = getName(obj)
           name=obj.name; 
        end
    end
    
end

