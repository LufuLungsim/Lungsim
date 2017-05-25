%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to write ASCII data files analogous to wBreath
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.2, 30. Januar 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok] = writeLogFileWbreath(input,output,table,parameters)

    minVentState	=   parameters.Simulation.minVentState;     % writes minute ventilations (per breath) into  logfile
    nBreathMax      =   parameters.Simulation.nBreathMax;       % maximal number of recorded breaths
    TOevaluation    =   parameters.Simulation.TOevaluation;
    lciStandard    	=   parameters.Simulation.lciStandard;
    lciByFit      	=   parameters.Simulation.lciByFit;
    tidal           =   parameters.Simulation.tidalAnalysis;       % tidal state
    tidalMean       =   parameters.Simulation.tidalMeanState;   % mean tidal data
    
    logFileName = 'logFile.txt';
    
    fid = fopen(logFileName, 'a');
    if fid<1
        error('writeLogFileWbreath:cannotOpenFile','cannot open file %s',logFileName);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Writing header
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ftell(fid)==0
        fprintf(fid,'Signal processing version\t');
        fprintf(fid,'data+time\tinput file name\toutput file name\tcalibration file name\tsettings file name\t');

        if tidal
            fprintf(fid,'minuteVentilationManeuvre [ml/min]\t');
            fprintf(fid,'breathRateManeuvre [#/min]\t');
            fprintf(fid,'Total number of breaths [-]\t');
        else
            fprintf(fid,'DeeBaseFRC [g/mol]\t');
            fprintf(fid,'DeeBaseLCI [g/mol]\t');
        end
        
        if  tidalMean
            if table.tidalMeans.breath0 > 0
                modifier = sprintf('%d-%d',table.tidalMeans.startBreath,table.tidalMeans.endBreath);
            else
                modifier = sprintf('(%d)-(%d)',table.tidalMeans.breath0-table.tidalMeans.nBreath+1,table.tidalMeans.breath0);
            end
            fprintf(fid,'MinuteVentilation%s [ml/min]\t',modifier);
            fprintf(fid,'BreathRate%s [1/min]\t',modifier);
            fprintf(fid,'inspDuration%s [s]\t',modifier);
            fprintf(fid,'expDuration%s [s]\t',modifier);
            fprintf(fid,'timeInspPeakFlow%s [s]\t',modifier);
            fprintf(fid,'timeExpPeakFlow%s [s]\t',modifier);
            fprintf(fid,'VolTidal%s [ml]\t',modifier);
            fprintf(fid,'tPTEF2tERatio%s [-]\t',modifier);
        end
        
        if lciStandard && ~tidal
            for i = 1:length(table.criticalEndRatios)
                fprintf(fid,'FRC(%1.1f%%) [l]\tLCI(%1.1f%%) [-]\tm0(%1.1f%%) [-]\tr1(%1.1f%%) [-]\tr2(%1.1f%%) [-]\t#Breaths(%1.1f%%) [-]\tduration(%1.1f%%) [s]\tminuteVentilation(%1.1f%%) [ml/min]\tbreathRate(%1.1f%%) [ml/min]\t', ...
                    table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100, ...
                    table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100);
            end
        end
        if TOevaluation && ~tidal
            for i = 1:length(table.criticalTurnovers)
                fprintf(fid,'FRCfromTO(%g) [l]\tDee(norm)fromTO(%g) [%%]\tm0fromTO(%g) [-]\tr1fromTO(%g) [-]\tr2fromTO(%g) [-]\t#BreathsForTO(%g) [-]\tdurationForTO(%g) [s]\t', ...
                    table.criticalTurnovers(i),table.criticalTurnovers(i),table.criticalTurnovers(i),table.criticalTurnovers(i),table.criticalTurnovers(i),...
                    table.criticalTurnovers(i),table.criticalTurnovers(i));
            end
        end
        if lciByFit && ~tidal
            for i = 1:length(table.criticalEndRatios)
                fprintf(fid,'FRCfit(%1.1f%%) [l]\tLCIfit(%1.1f%%) [-]\tm0Fit(%g) [-]\tr1Fit(%g) [-]\tr2Fit(%g) [-]\t#BreathsFit(%1.1f%%) [-]\tdurationFit(%1.1f%%) [s]\tminuteVentilationFit(%1.1f%%) [ml/min]\tbreathRateFit(%1.1f%%) [ml/min]\t', ...
                    table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,...
                    table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100,table.criticalEndRatios(i)*100);
            end
        end
        fprintf(fid,'VolInspMean [ml]\tVolExpMean [ml]\trelVolInspStd [%%]\trelVolExpStd [%%]\t');
        
        if minVentState && ~tidal
            for i = 1:nBreathMax
                fprintf(fid,'inspiredVolume(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'expiredVolume(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'meanInspiredFlow(%d) [ml/min]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'meanExpiredFlow(%d) [ml/min]\t',i);
            end
        end
        
        if minVentState && tidal
            for i = 1:nBreathMax
                fprintf(fid,'Inspiration time(%d) [s]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Expiration time(%d) [s]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Total breath time(%d) [s]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Peak tidal inspiratory flow(%d) [ml/s]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Peak tidal expiratory flow(%d) [ml/s]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Time to PTIF(%d) [s]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Time to PTEF(%d) [s]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Inspired Volume(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Expired Volume(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Tidal Volume(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Respiratory Rate(%d) [1/min]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Ratio tPTEF/tE(%d) [%%]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Mean Inspired Flow(%d) [ml/min]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Mean Expired Flow(%d) [ml/min]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Minute Ventilation(%d) [ml/min]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'Breath Typ(%d) [wi|wo|pre]\t',i);
            end
        end
        
        fprintf(fid,'\n');
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Writing data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nBreaths=floor((table.nHalfBreaths+1)/2);
    
    calibration = parameters.Calibration;
    simulation  = parameters.Simulation;
    
    [pathString,nameInput,extension]=fileparts(input);
    nameInput=strcat(nameInput,extension);
    [pathString,nameOutput,extension]=fileparts(output);
    nameOutput=strcat(nameOutput,extension);
    [pathString,nameCalibration,extension]=fileparts(calibration.filePath);
    nameCalibration=strcat(nameCalibration,extension);
    [pathString,nameSettings,extension]=fileparts(calibration.settingsFilePath);
    nameSettings=strcat(nameSettings,extension);
    
    fprintf(fid,'%s\t',simulation.versionString);
    fprintf(fid,'%s\t%s\t%s\t%s\t%s\t',datestr(now),nameInput,nameOutput,nameCalibration,nameSettings);
    
    if tidal
        fprintf(fid,'%e\t',table.minuteVentilationTotal(1)*1e6*60);
        fprintf(fid,'%e\t',table.breathRateTotal(1)*60);
        % % % % % fprintf(fid,'%e\t',table.minuteVentilationCrit(1)*1e6*60);
        % % % % % fprintf(fid,'%e\t',table.breathRateCrit(1)*60);
        fprintf(fid,'%d\t',nBreaths);
    else
        fprintf(fid,'%e\t',table.DeeBaseFRC*1000);
        fprintf(fid,'%e\t',table.DeeBaseLCI*1000);
    end
    
    if tidalMean
        ref=table.tidalMeans;
        fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t',ref.minuteVentilation*1e6*60,ref.breathRate*60,ref.inspDuration,ref.expDuration,ref.timeInspPeakFlow,ref.timeExpPeakFlow,ref.VolTidal*1e6,ref.tPTEF2tERatio);
    end
    
    if lciStandard && ~tidal
        Data=table.Direct;
        for i = 1:length(table.criticalEndRatios)
            fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%g\t%e\t%e\t%e\t',Data.FRC(i)*1e3,Data.LCI(i),Data.m0(i),Data.r1(i),Data.r2(i),Data.nBreaths(i),Data.duration(i),Data.minuteVentilation(i)*1e6*60,Data.breathRate(i)*60);
        end
    end
    if TOevaluation && ~tidal
        Data=table.Turnover;
        for i = 1:length(table.criticalTurnovers)
            fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%g\t%e\t',        Data.FRC(i)*1e3,Data.DeeNorm(i),Data.m0(i),Data.r1(i),Data.r2(i),Data.nBreaths(i),Data.duration(i));
        end
    end
    if lciByFit && ~tidal
        Data=table.Fit;
        for i = 1:length(table.criticalEndRatios)
            fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%g\t%e\t%e\t%e\t',Data.FRC(i)*1e3,Data.LCI(i),Data.m0(i),Data.r1(i),Data.r2(i),Data.nBreaths(i),Data.duration(i),Data.minuteVentilation(i)*1e6*60,Data.breathRate(i)*60);
        end
    end
    
    fprintf(fid,'%e\t%e\t',table.VolInspMean*1e6,table.VolExpMean*1e6,table.VolInspStd*100,table.VolExpStd*1e2);
    
    if minVentState && ~tidal
        nBreathEff=min(nBreaths,nBreathMax);
        nRest=max(0,nBreathMax-nBreathEff);
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.VolInsp(j)*1e6);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.VolExp(j)*1e6);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.FlowInspMean(j)*1e6*60);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.FlowExpMean(j)*1e6*60);
        end
        writeTabs(fid,nRest)
    end
    
    if minVentState && tidal
        nBreathEff=min(nBreaths,nBreathMax);
        nRest=max(0,nBreathMax-nBreathEff);
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.inspDuration(j));
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.expDuration(j));
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.breathDuration(j));
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.inspPeakFlow(j)*1e6);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.expPeakFlow(j)*1e6);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.timeInspPeakFlow(j));
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.timeExpPeakFlow(j));
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.VolInsp(j)*1e6);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.VolExp(j)*1e6);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.VolTidal(j)*1e6);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',60/table.breathDuration(j));
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.tPTEF2tERatio(j)*100);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.FlowInspMean(j)*1e6*60);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',table.FlowExpMean(j)*1e6*60);
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%e\t',60*table.VolTidal(j)*1e6/table.breathDuration(j));
        end
        writeTabs(fid,nRest)
        
        for j = 1:nBreathEff
            fprintf(fid,'%s\t',table.breathType{j});
        end
        writeTabs(fid,nRest)
    end
    
    fprintf(fid,'\n');
    fclose(fid);
    
    ok=1;
end % writeLogFileWbreath


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find index in list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function []=writeTabs(fid,n)

    for k = 1:n
        fprintf(fid,'\t');
    end
end

