function [  ] = savePredictionData(parameters, table, gas, outputModel,  outputMeasurement, outputModelInterp, outputMeasurementInterp, LCIStopCrit)

    fileInput = parameters.Simulation.fileInput;
    [path, name, extension] = fileparts(fileInput);
       
    fileID = ['..\' name '.xlsx'];
    
    while xlswrite(fileID, 'test') == 0
        waitfor(msgbox({'Can not continue because following file is open';  fileID; '' ; 'Please close the file.'}));
    end
    
    header = {'LCI',	'Breath 2_5',	'CEV',	'VN2Exp',	'FRC'};
    values = {gas.Standard.lci(1), gas.Standard.nBreaths(1), gas.General.cev_ds(gas.Standard.nBreaths(1))*1000, gas.General.cumNetExpVol(gas.Standard.nBreaths(1))*1000, gas.General.frcao(gas.Standard.nBreaths(1))*1000};
    
    xlswrite(fileID, [header; values], 'Base data', 'A1');
    xlswrite(fileID, LCIStopCrit, 'Base data', 'B4');
    xlswrite(fileID, {'LCIasDefinedByStopCriteria'; 'Rel error of prediction'; 'Stop Criteria Breath'; 'Stop Criteria Breath Ratio to LCI breath'}, 'Base data', 'A6');
    
    %xlswrite(fileID,{'Values interpolated based purely on functions'}, 'Interpolated data', 'A1');
    %xlswrite(fileID,outputModelInterp, 'Interpolated data', 'A2');
    
    offset = size(outputModelInterp, 1);
    xlswrite(fileID,{'Values interpolated based on functions and measurement values'}, 'Interpolated data', ['A1']);
    xlswrite(fileID,outputMeasurementInterp, 'Interpolated data', ['A2']);

    %xlswrite(fileID,{'Values Calculated based purely on functions'}, 'Breath data', 'A1');
    %xlswrite(fileID,outputModel, 'Breath data', 'A2');
    
    offset = size(outputModel, 1);
    xlswrite(fileID,{'Values Calculated based on functions and measurement values'}, 'Breath data', ['A1']);
    xlswrite(fileID,outputMeasurement, 'Breath data', ['A2']);
    
end

