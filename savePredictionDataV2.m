function [  ] = savePredictionData(parameters, table, gas, outputModel, outputModelInterp, outputParams, LCIStopCrit)
 

    baseDataAnnotation = {'LCIasDefinedByStopCriteria'; 'Rel error of prediction'; 'CEV'; 'FRC'; 'M0'; 'M1'; 'M2';'Stop Criteria Breath'; 'Stop Criteria Breath Ratio to LCI breath'};
    
    if isfield(gas.Standard, 'momentRatio') && length(gas.Standard.momentRatio) > 1
        header = {'LCI',	'Breath 2_5',	'CEV',	'VN2Exp',	'FRC', 'M0', 'M1', 'M2'};
        values = {gas.Standard.lci(1)                                       , ...
                  gas.Standard.nBreaths(1)                                  , ...
                  gas.General.cev_ds(gas.Standard.nBreaths(1))*1000         , ...
                  gas.General.cumNetExpVol(gas.Standard.nBreaths(1))*1000   , ...
                  gas.General.frcao(gas.Standard.nBreaths(1))*1000          , ...
                  gas.Standard.momentRatio(1,1)                             , ...
                  gas.Standard.momentRatio(1,2)                             , ...
                  gas.Standard.momentRatio(1,3)};
    else 
        header = {'LCI',	'Breath 2_5',	'CEV',	'VN2Exp',	'FRC'};
        values = {gas.Standard.lci(1)                                       , ...
                  gas.Standard.nBreaths(1)                                  , ...
                  gas.General.cev_ds(gas.Standard.nBreaths(1))*1000         , ...
                  gas.General.cumNetExpVol(gas.Standard.nBreaths(1))*1000   , ...
                  gas.General.frcao(gas.Standard.nBreaths(1))*1000};
    end
    
    fileInput = parameters.Simulation.fileInput;
    [path, name, extension] = fileparts(fileInput);
       
    fileID = [pwd '\SimulationData\' name '.xlsx'];
    %Connect to excel
    e = actxserver('Excel.Application');
    % Get Workbook object
    eWorkbook = e.Workbooks.Add;
    %e.Visible = 1;
    %Get Worksheets object
    eSheets = eWorkbook.Sheets;
    % Add after the last sheet
    eSheets.Add([], eSheets.Item(eSheets.Count));
    %Get activate and name the sheet
    eSheet1 = eSheets.get('Item',2);
    eSheet1.Activate
    eSheet1.Name = 'Base data';
    
    %Save data to Base Data
    eActivesheetRange = get(e.Activesheet,'Range',['A1:' ExcelCol(size(header,2)) '2']);
    eActivesheetRange.Value = [header; values];
    
    eActivesheetRange = get(e.Activesheet,'Range',['B4:' ExcelCol(size(LCIStopCrit,2)+1) '' num2str(size(LCIStopCrit,1)+3)]);
    eActivesheetRange.Value = LCIStopCrit;
    
    eActivesheetRange = get(e.Activesheet,'Range',['A6:A' num2str(5+length(baseDataAnnotation))]);
    eActivesheetRange.Value = baseDataAnnotation;

    % Add after the last sheet
    eSheets.Add([], eSheets.Item(eSheets.Count));
    %Get activate and name the sheet
    eSheet2 = eSheets.get('Item',3);
    eSheet2.Activate
    eSheet2.Name = 'Interpolated Data';
    
    %Save data
    eActivesheetRange = get(e.Activesheet,'Range',['A1']);
    eActivesheetRange.Value = {'Values interpolated based on functions and measurement values'};
     
    eActivesheetRange = get(e.Activesheet,'Range',['A2:' ExcelCol(size(outputModelInterp,2)) '' num2str(size(outputModelInterp,1)+1)]);
    eActivesheetRange.Value = outputModelInterp;
    
    % Add after the last sheet
    eSheets.Add([], eSheets.Item(eSheets.Count));
    %Get activate and name the sheet
    eSheet3 = eSheets.get('Item',4);
    eSheet3.Activate
    eSheet3.Name = 'Breathwise Data';
    
    %Save data
    eActivesheetRange = get(e.Activesheet,'Range',['A1']);
    eActivesheetRange.Value = {'Values Calculated based on functions and measurement values'};
     
    eActivesheetRange = get(e.Activesheet,'Range',['A2:' ExcelCol(size(outputModel,2)) '' num2str(size(outputModel,1)+1)]);
    eActivesheetRange.Value = outputModel;
    
    % Add after the last sheet
    eSheets.Add([], eSheets.Item(eSheets.Count));
    %Get activate and name the sheet
    eSheet4 = eSheets.get('Item',5);
    eSheet4.Activate
    eSheet4.Name = 'Prediction Params';
    
    %Save data
    eActivesheetRange = get(e.Activesheet,'Range',['A1']);
    eActivesheetRange.Value = {'Parameters Used for Prediction'};
     
    eActivesheetRange = get(e.Activesheet,'Range',['A2:' ExcelCol(size(outputParams,2)) '' num2str(size(outputParams,1)+1)]);
    eActivesheetRange.Value = outputParams;
        
    %Save and close the datasheet
    SaveAs(eWorkbook,fileID)
    eWorkbook.Saved = 1;
    Close(eWorkbook)
    
    %Quit and close the spreadsheet
    Quit(e)
    delete(e)

%     while xlswrite(fileID, 'test') == 0
%         waitfor(msgbox({'Can not continue because following file is open';  fileID; '' ; 'Please close the file.'}));
%     end
%     
%     
%     xlswrite(fileID, [header; values], 'Base data', 'A1');
%     xlswrite(fileID, LCIStopCrit, 'Base data', 'B4');
%     xlswrite(fileID, {'LCIasDefinedByStopCriteria'; 'Rel error of prediction'; 'CEV'; 'FRC'; 'M0'; 'M1'; 'M2';'Stop Criteria Breath'; 'Stop Criteria Breath Ratio to LCI breath'}, 'Base data', 'A6');
%     
%     xlswrite(fileID,{'Values interpolated based on functions and measurement values'}, 'Interpolated data', ['A1']);
%     xlswrite(fileID,outputModelInterp, 'Interpolated data', ['A2']);
% 
%     xlswrite(fileID,{'Values Calculated based on functions and measurement values'}, 'Breath data', ['A1']);
%     xlswrite(fileID,outputModel, 'Breath data', ['A2']);
%     
%     xlswrite(fileID,{'Parameters Used for Prediction'}, 'Prediction Params', ['A1']);
%     xlswrite(fileID,outputParams, 'Prediction Params', ['A2']);   
end

function  Out=ExcelCol(n)

    %Optional to change representation and base
    ABC='ABCDEFGHIJKLMNOPQRSTUVWXYZ';
    base=26; 
    diff=1;
    i=0;
    
    while diff<=n
        letter_ind=1+mod(floor((n-diff)/base^i),base);
        i=i+1;
        temp(i)=ABC(letter_ind);
        diff=diff+base^i;
    end

    Out = temp;
end

