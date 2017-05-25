function [ table ] = getGasFromSignal( signal )
    
    %A List of possible gases which could be contained in the signal vector
    possibleGases={'N2', 'O2', 'CO2', 'Ar', 'He', 'SF6', 'MMss', 'MMms'};
    for i=1:length(possibleGases)
        %If the gas is in the signal vector, create a class instance which
        %is saved in the table struct under the gas name and can be
        %identified by its name.         
        if isfield(signal,possibleGases{i});
            table.(matlab.lang.makeValidName(possibleGases{i}))=Gas();
            table.(matlab.lang.makeValidName(possibleGases{i})).name=possibleGases{i};
        end
    end
end

