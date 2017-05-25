%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Apply signal correction to breathing maneuver data files
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.1, 28. August 2012
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters]=applyCorrectionSpirowareMassSpec(sourcePathName,fileName,destinationPathName,parameters,flag)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set general parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    verb                =   parameters.Simulation.verb;                     % verbosity state
    
    bfileOut            =   parameters.Simulation.bfileOutput;              % bfile (B-) output
    cfileOut            =   parameters.Simulation.cfileOutput;              % cfile (C-) output
    inselOut            =   parameters.Simulation.inselOutput;              % insel (I-) output
    lfileOut            =   parameters.Simulation.leicesterOutput;          % leicester (L-) output
    
    perFileCalibration  =   parameters.Simulation.fileCalibrationYes;       % per file calibration is set
    modifier            =   parameters.Simulation.fileCalibrationModifier;  % per file modifier
    
    manualCalibration   =   parameters.Simulation.manualState;              % manual calibration
    autoCalibration     =   parameters.Simulation.autoState;                % automatic calibration
    cropping            =   parameters.Simulation.croppingState;            % data cropping
    recalibrationState  =   parameters.Operation.recalibration;             % recalibration state
    torontoFile         =   parameters.Simulation.torontoFile;              % toronto data type
    onlyPlotState       =   parameters.Simulation.onlyPlotState;            % only raw signals are plotted
    aFile               =   parameters.Simulation.aFile;
    bFile               =   parameters.Simulation.bFile;
    
    if torontoFile
        parameters.Simulation.dt = 0.01;                                    % has different time scale
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filename as string
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~iscell(fileName)
        fileName={fileName};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform correction and write various file formats to disk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    outputType=cfileOut+inselOut*2+bfileOut*4+lfileOut*8+torontoFile*16;
    outputChar = getOutputChar(outputType);
  
    hWaitbar = findall(0,'Tag','TMWWaitbar');
    if isempty(hWaitbar)
        hWaitbar = waitbar(0,'Starting Process', 'Name' ,'Processing Files'); 
    end
    numOfFile = length(fileName);
    for i=1:numOfFile
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         % Calibration
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         if ~bFile
%             if autoCalibration
%                 parameters = fileCalibration(sourcePathName, i, fileName, recalibrationState, parameters, verb, 2);
%             elseif manualCalibration
%                 parameters = fileCalibration(sourcePathName, i, fileName, recalibrationState, parameters, verb, 1);
%             elseif perFileCalibration
%                 parameters = calibrationPerFile(i, fileName, sourcePathName, postfix, divider,recalibrationState, parameters);
%             end 
%         end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Load signal Data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fileNameString = cell2mat(fileName(i));
        parameters = setSpirowareType(fileNameString,parameters);
        fileInput=fullfile(sourcePathName,fileNameString);
        parameters.Simulation.fileInput = fileInput;
        fprintf('\nsignal correction for file: %s\n',fileInput)
        [header, data, signal]=readSignalSpirowareMassSpec(fileInput,i,parameters);
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Crop and Plot data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if onlyPlotState
            plotSignalSpirowareMassSpec(signal,parameters,i*flag);
            continue
        end
        if cropping
            [time.min,time.max]=cropDataSpirowareMassSpec(signal,parameters,'data cropping');
            signal = cropSignal(time, signal);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Signal correction 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [signalCorrected,breathIndex]   =   signalCorrection(signal,parameters,i*flag);
        waitbar((i-1)/numOfFile+1/numOfFile*0.4, hWaitbar, 'Signal Corrected'); 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Segment signal to breaths
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        table   =   getBreathTableMBW(breathIndex,signalCorrected,parameters);   % breath table generation
        if ~torontoFile
            meanN2  =   mean(signalCorrected.N2Poi3);
            stdN2   =   std(signalCorrected.N2Poi3);
            fprintf('meanN2=%g, stdN2=%g, errorN2=%g %%\n' ,meanN2,stdN2,stdN2/meanN2*100)
        end
        
        waitbar((i-1)/numOfFile+1/numOfFile*0.6, hWaitbar, 'Signal Analyzed');
        
        [nameSignal, nameTable] = nameSignalAndTable(fileName, i, outputChar);
                
        fileOutput=fullfile(destinationPathName,nameTable);
        fprintf('\nwriting table to file: %s\n',fileOutput)
        
        %Write Breath table (Get values from signal and write to a table (Calculation of lci & FRC)
        [table, speciesNames]   = writeBreathTableMBWWrapper(table,fileOutput,parameters,flag*i, signalCorrected);
        
        waitbar((i-1)/numOfFile+1/numOfFile*0.8, hWaitbar, 'Fields Calculated');
        
        if lfileOut % anonymize file name
            [nameSignal,parameters.Simulation.nameHash]=nameAnonymizer(nameSignal,parameters.Simulation.nameHash);
        end
        
        fileOutput=fullfile(destinationPathName,nameSignal);
        fprintf('\nwriting corrected signal to file: %s\n\n',fileOutput)
        writeSignalSpirowareMassSpec(signalCorrected,fileOutput,outputType,verb);
        
        for n = 1:length(speciesNames)
            writeLogFileMBW(fileInput,fileOutput,table,speciesNames(n),parameters); 
        end  
    end
    
    close(hWaitbar)
end

function outputChar = getOutputChar(outputType)
    switch outputType
        case 1
            outputChar = 'C';
        case 2
            outputChar = 'I';
        case 4
            outputChar = 'B';
        case 8
            outputChar = 'L';
        case 16
            outputChar = 'T';
        otherwise
            outputChar = '';
    end
end

function [table, speciesNames]   = writeBreathTableMBWWrapper(table,fileOutput,parameters,flag, signalCorrected)
        
    torontoFile        =   parameters.Simulation.torontoFile;
    dtgmbw             =   parameters.Simulation.dtgmbwAnalysis;
    dtgmbwBaby         =   parameters.Simulation.dtgmbwBaby;     %SF6 4% %%TODO rename dtgmbwBaby to dtgmbw4SF6
    bFile              =   parameters.Simulation.bFile;
    
    speciesNames = [];
    
    if torontoFile
        speciesNames{1}='He';
        speciesNames{2}='SF6';
        
        table	= writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames, 0, 0);
    else
        
        %SF6 4%
        if dtgmbwBaby
            if bFile
                speciesNames{1}='SF6';
                if isWashin(table.MMss.et) %O2 Washin
                    table   = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 1, 1); 
                    speciesNames{1}='SF6WI';    
                    table.SF6WI = table.SF6;
                else %O2 Washout
                    table   = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 0, 0);                                
                    speciesNames{1}='SF6WO';    
                    table.SF6WO = table.SF6;
                end
            else % 4% SF6 Washin Washout for A file
                speciesNames{1}='SF6';
                tableSF6WI = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 1, 1);
                tableSF6WO = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 0, 0);
                                
                speciesNames{1}='SF6WI';      
                speciesNames{2}='SF6WO';
                
                table = tableSF6WI;
                table.SF6WI = tableSF6WO.SF6;
                table.SF6WO = tableSF6WI.SF6;
            end              
        %SF6 3% O2 97%
        elseif dtgmbw
            speciesNames{1}='N2';   
            speciesNames{2}='SF6';
            if bFile
                if isWashin(table.O2.et) %O2 Washin
                    tableSF6WI   = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(2), 1, 1);
                    tableN2WO    = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 1, 0);
                    
                    speciesNames{1}='N2WO';
                    speciesNames{2}='SF6WI';
                    table = tableN2WO;
                    table.N2WO  = tableN2WO.N2;
                    table.SF6WI = tableSF6WI.SF6; 
                else %O2 Washout
                    tableSF6WO   = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(2), 0, 0);
                    tableN2WI    = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 0, 1);
                    
                    speciesNames{1}='N2WI';
                    speciesNames{2}='SF6WO';
                    table = tableN2WI;
                    table.N2WI  = tableN2WI.N2;
                    table.SF6WO = tableSF6WO.SF6;                         
                end
            else % 3% SF6 Washin Washout for A file
                tableSF6WI = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(2), 1, 1);
                tableSF6WO = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(2), 0, 0);
                tableN2WI  = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 0, 1);
                tableN2WO  = writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 1, 0);                
                
                speciesNames{1}='N2WO';
                speciesNames{2}='SF6WI';              
                speciesNames{3}='N2WI';
                speciesNames{4}='SF6WO';
                
                table = tableN2WI;
                table.N2WO  = tableN2WO.N2;
                table.SF6WI = tableSF6WI.SF6;
                table.N2WI  = tableN2WI.N2;
                table.SF6WO = tableSF6WO.SF6;
            end           
        else
            speciesNames{1}='N2';   
            table	= writeBreathTableMBW(table,fileOutput,parameters,flag, signalCorrected, speciesNames(1), 1, 0);
        end

    end
end

function p = isWashin(concentration)
    p = polyfit(1:length(concentration), concentration', 1);
    p = p(1)>0;
end

function [divider, postfix] = extractPossibleCalibrationFileNameQualifiers(modifier)
    indexFind=findstr(';',modifier);
    if ~isempty(indexFind)
        if length(indexFind)>1
            endString=indexFind(2)-1;
        else
            endString=length(modifier);
        end
        divider=modifier(1:indexFind(1)-1);
        postfix=modifier(indexFind(1)+1:endString);
    else
        divider='';
        postfix=modifier;
    end
end

function calibrationFileName=getCalibrationFileName(i, fileName, sourcePathName, postfix, divider)
    calibrationFileName=cell2mat(fileName(i));
    [pathString,stem,extension]=fileparts(calibrationFileName);
    pattern=strcat(sourcePathName,stem,postfix,extension);
    listing=ls(pattern);
    if ~isempty(listing)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % first pass: try to match full file name (with postfix)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        calibrationFileName=strtrim(listing(1,:));
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % second pass: try to match file name with alternative variant number (with postfix)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        indexMatch=findstr(divider,calibrationFileName);
        if isempty(indexMatch)
            warning('(applyCorrectionSpirowareMassSpec): file name does not contain pattern to associate calibration file! using existing parameter instead');
            calibrationFileName='';
        else
            indexMatch=indexMatch(max(1,end-1));
            stem=calibrationFileName(1:indexMatch-1);
            pattern=strcat(sourcePathName,stem,divider,'*',postfix,extension);
            listing=ls(pattern);
            if isempty(listing)
                warning('(applyCorrectionSpirowareMassSpec): no calibration file with pattern found! using existing parameter instead');
                calibrationFileName='';
            else
                calibrationFileName=strtrim(listing(1,:));
            end
        end
    end
end

function signal = cropSignal(time, signal)
    imin    =   find(signal.ts>=time.min,1);
    imax    =   find(signal.ts<=time.max,1,'last');

    fields  =   fieldnames(signal);
    for j=1:length(fields)
        if length(signal.(fields{j}))>=imax
            signal.(fields{j})  =   signal.(fields{j})(imin:imax);
        end
    end
end

function [nameSignal, nameTable] = nameSignalAndTable(fileName, i, outputChar)
    nameSignal=cell2mat(fileName(i));
    nameTable =cell2mat(fileName(i));
    if strcmp('A-',nameSignal([1,2]))
        nameSignal(1)=outputChar;
        nameTable(1)='e';
        nameTable=strcat('Tabl',nameTable);
    else
        nameSignal=strcat(outputChar,'-',nameSignal);
        nameTable =strcat('Table-',nameTable);
    end
end

function parameters = calibrationPerFile(i, fileName, sourcePathName, postfix, divider,recalibrationState, parameters)
    calibrationFileName = getCalibrationFileName(i, fileName, sourcePathName, postfix, divider);           
    if ~isempty(calibrationFileName)
        filePath=fullfile(sourcePathName,calibrationFileName);
        [parameters]=getCorrectionCoeff(recalibrationState,filePath,parameters,-1,0);
    end
end

function parameters = fileCalibration(sourcePathName, i, fileName, recalibrationState, parameters, verb, type)
    filePath=fullfile(sourcePathName,cell2mat(fileName(i)));
    parameters = getCorrectionCoeff(recalibrationState,filePath,parameters,type,verb);
end
