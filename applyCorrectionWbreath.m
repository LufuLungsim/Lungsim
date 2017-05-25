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
function [parameters]=applyCorrectionWbreath(sourcePathName,fileName,destinationPathName,parameters,flag)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set general parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cropping            =   parameters.Simulation.croppingState;            % data cropping
    onlyPlotState       =   parameters.Simulation.onlyPlotState;            % only raw signals are plotted
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % filename as string
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ~iscell(fileName)
        fileName={fileName};
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Perform correction and write various file formats to disk
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=1:length(fileName)
        fileNameString = cell2mat(fileName(i));
        fileInput=fullfile(sourcePathName,fileNameString);
        fprintf('\nsignal correction for file: %s\n',fileInput)
        [header data signal]=readSignalWbreath(fileInput,i,parameters);
        
        if cropping
            [time.min,time.max]=cropDataWbreath(signal,parameters,'data cropping');
            imin    =   find(signal.ts>=time.min,1);
            imax    =   find(signal.ts<=time.max,1,'last');
            
            fields  =   fieldnames(signal);
            for j=1:length(fields)
                signal.(fields{j})  =   signal.(fields{j})(imin:imax);
            end
        end
        
        if onlyPlotState
            plotSignalWbreath(signal,parameters,i*flag);
            continue
        end
        
        [signalCorrected,breathIndex]=signalCorrectionWbreath(signal,parameters,i*flag);
        
        table=getBreathTableWbreath(breathIndex,signalCorrected,parameters);
        
        fileOutput=fullfile(destinationPathName,fileNameString);
        fprintf('\nwriting table to file: %s\n',fileOutput)
        [ok,table]=writeBreathTableWbreath(table,fileOutput,parameters,flag*i);
        
        writeLogFileWbreath(fileInput,fileOutput,table,parameters);
        
    end
end
