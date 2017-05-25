%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Get spiroware data type of A-files
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.1, 17. Sept. 2013
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters]=setSpirowareType(sourceFileName,parameters)

    findArray = strfind(sourceFileName,'A-');
    findBArray= strfind(sourceFileName,'B-');
    if isempty(findArray)&&isempty(findBArray)
        parameters.Operation.spirowareType = -1;
        infoString = 'data file not of Spiroware type';
    elseif ~isempty(findBArray)
        infoString = 'B-File';
        parameters.Operation.spirowareType = 3;
    else
        startDate = findArray(1)+2;
        dateString=sourceFileName(startDate:startDate+7);
        if  isnan(str2double(dateString))
        	parameters.Operation.spirowareType = -1;
            infoString = sprintf('cannot determine date from string ''%s''',dateString);
        else
            dateNumber=datenum(str2double(dateString(1:4)),str2double(dateString(5:6)),str2double(dateString(7:8)));
            if dateNumber < parameters.Operation.dateCO2CORR
                parameters.Operation.spirowareType = 2;
                infoString = 'CO2 signal corrected to dry conditions';
            elseif dateNumber >= parameters.Operation.dateCO2RAW
                parameters.Operation.spirowareType = 0;
                infoString = 'raw CO2 signal';
            else
                parameters.Operation.spirowareType = 1;
                infoString = 'raw CO2 signal, but not named _RAW';
            end
        end
    end
    
    fprintf('Sprioware data file type = %d (%s)\n', parameters.Operation.spirowareType, infoString);
end
