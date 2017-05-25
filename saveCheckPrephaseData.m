function [] = saveCheckPrephaseData( table, range, parameters )

    fileName = parameters.Simulation.fileInput;

    out = range';
    out = [out, table.N2.et(range)*100];
    out = [out, table.N2.eti(range)*100];
    out = [out table.N2.exp(range)*1000];
    out = [out table.N2.insp(range)*1000];
    out = [out table.O2.et(range)*100];
    out = [out table.O2.eti(range)*100];
    out = [out table.O2.exp(range)*1000];
    out = [out table.O2.insp(range)*1000];
    out = [out table.CO2.et(range)*100];
    out = [out table.CO2.eti(range)*100];
    out = [out table.CO2.exp(range)*1000];
    out = [out table.CO2.insp(range)*1000];
    out = [out table.volExp(range)*1000];
    out = [out table.volInsp(range)*1000];
    out = out';
    out = mat2cell(out, ones(size(out, 1),1), ones(size(out,2),1));
    
    rowNames = {'Breath Number'; ...
                'N2 end expiratory concentration'; 'N2 end inspiratory concentration'; 'N2 expiratory volume'; 'N2 inspiratory volume'; ...
                'O2 end expiratory concentration'; 'O2 end inspiratory concentration'; 'O2 expiratory volume'; 'O2 inspiratory volume'; ...
                'CO2 end expiratory concentration';'CO2 end inspiratory concentration';'CO2 expiratory volume';'CO2 inspiratory volume'; ...
                'Expiratory volume'; 'Inspiratory volume'};
            
    out = [rowNames, out];
    
    folder   = pwd;
    saveFile = 'PrephaseInfo.txt';
    fileID   = fopen([folder '\' saveFile], 'a');
    fprintf(fileID, '%s \r\n', fileName);
    
    formatSpec = ['%s; ' repmat('%0.3f; ', 1, length(range)) '\r\n'];
    for row = 1:length(rowNames)
    	fprintf(fileID, formatSpec, out{row,:});        
    end
    
    fprintf(fileID, '\r\n');
    fclose(fileID);
end

