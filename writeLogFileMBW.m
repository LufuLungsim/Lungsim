%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to write ASCII data files analogous to wBreath
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.1, 28. August 2012
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok] = writeLogFileMBW(input,output,table,speciesNames,parameters)

    slopesState     =   parameters.Simulation.slopesState;      % writes slopes (per breath) into  logfile
    minVentState	=   parameters.Simulation.minVentState;     % writes minute ventilations (per breath) into  logfile
    tidal           =   parameters.Simulation.tidalAnalysis;    % contains data from tidal analysis
    tidalMean       =   parameters.Simulation.tidalMeanState;   % contains mean tidal data
    n2mbw           =   parameters.Simulation.n2mbwAnalysis;    % contains data from N2 MBW analysis
    dtgsbw          =   parameters.Simulation.dtgsbwAnalysis;   % contains data from DTG-SBW analysis
    nBreathMax      =   parameters.Simulation.nBreathMax;       % maximal number of recorded breaths
    ScondSacin      =   parameters.Simulation.ScondSacin;       % evaluation of Scond/Sacin
    ScondSacinNorm  =   parameters.Simulation.ScondSacinNorm;   % Normalization type of Scond/Sacin 
    capnoIndices    =   parameters.Simulation.CapnoIndices;
    TOevaluation    =   parameters.Simulation.TOevaluation;
    lciStandard    	=   parameters.Simulation.lciStandard;
    lciByFit      	=   parameters.Simulation.lciByFit;
    dtgmbw          =   parameters.Simulation.dtgmbwAnalysis;
    torontoFile     =   parameters.Simulation.torontoFile;
    momentRatio     =   parameters.Simulation.MomentRatio;
    fastSlow        =   parameters.Simulation.fastSlow;
    
    nGases          =   length(speciesNames);
    
    if parameters.Simulation.fullOutput
        
    if dtgsbw
        nBreathMax = 1; % to reduced the output column count (otherwise empty cells)
    end
    
    logFileName = 'logFile.txt';
    
    fid = fopen(logFileName, 'a');
    if fid<1
        error('cannot open file %s\n',logFileName);
    end
    
    if ftell(fid)==0
        fprintf(fid,'Signal processing version\t');
        if torontoFile
            fprintf(fid,'data+time\tinput file name\toutput file name\tsettings file name\t');
        else
            fprintf(fid,'data+time\tinput file name\toutput file name\t');
            if n2mbw
                fprintf(fid,'calibration file name\t');
            end
            fprintf(fid,'settings file name\t');
            if n2mbw
                fprintf(fid,'delay Iv->CO2\trecalibration\t');
                %fprintf(fid,'delay Iv->CO2\tdelay CO2->O2\tdelay CO2->MMss\ttauCO2\tkCO2\tkO2\tkMMss\trecalibration\t');
            end
        end
        fprintf(fid,'BTPS\tstatic\tdynamic\t');
        fprintf(fid,'BTPSinsp\tBTPSexp\t');
        if dtgmbw
            fprintf(fid,'Analysed Gas\t');
        end
        if dtgsbw
            fprintf(fid,'SlopeCO2 [%%/l]\tSlopeN2 [%%/l]\tSlopeMMss [g/mol/l]\tSlopeMMssCalc [g/mol/l]\tSlopeDiffMMss [g/mol/l]\t');
            fprintf(fid,'InterceptionVolume [%%/ExpVolume]\tEndTidalDiffMMss [g/mol]\tDiffMMssMax [g/mol]\tDiffMMssMaxVolume [%%ExpVol]\t');
            fprintf(fid,'ExpiredVolume [l]\t');
            fprintf(fid,'minuteVentilationManeuvre [ml/min]\t');
            fprintf(fid,'breathRateManeuvre [#/min]\t');
        end
        
        if n2mbw || dtgmbw
            for j = 1 : nGases
                if torontoFile
                    gasName=strcat('(',speciesNames{j},')');
                else
                    gasName='';
                end
                
                if ScondSacinNorm
                    normType = '';
                else
                    normType = '*VTexp';
                end
                
                if ScondSacin
                    if ScondSacinNorm
                        fprintf(fid,'Scond%s [1/l/TO]\tSacin%s [1/l]\tR2%s [-]\t#Breaths total%s [-]\t#Breaths excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName);
                        fprintf(fid,'ScondStar%s [1/l/TO]\tSacinStar%s [1/l]\tR2%s [-]\t#BreathsStar total%s [-]\t#BreathsStar excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName);
                    else
                        fprintf(fid,'Scond%s [1/TO]\tSacin%s [-]\tR2%s [-]\t#Breaths total%s [-]\t#Breaths excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName);
                        fprintf(fid,'ScondStar%s [1/TO]\tSacinStar%s [-]\tR2%s [-]\t#BreathsStar total%s [-]\t#BreathsStar excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName);
                    end
                end
                
                if torontoFile
                    gasName=strcat(speciesNames{j},':');
                else
                    gasName='';
                end
                if lciStandard
                    for i = 1:length(table.criticalEndRatios)
                        endRatio=table.criticalEndRatios(i)*100;
                        fprintf(fid,'FRC(%s%1.1f%%) [l]\tLCI(%s%1.1f%%) [-]\t#Breaths(%s%1.1f%%) [-]\tduration(%s%1.1f%%) [s]\tminuteVentilation(%s%1.1f%%) [ml/min]\tbreathRate(%s%1.1f%%) [1/min]\t', ...
                            gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio);
                        if momentRatio
                            fprintf(fid, 'mr0(%s%1.1f%%)\tmr1(%s%1.1f%%)\tmr2(%s%1.1f%%)\t', ...
                                         gasName,endRatio,gasName,endRatio,gasName,endRatio);
                        end
                    end
                end
                if lciByFit
                    for i = 1:length(table.criticalEndRatios)
                        endRatio=table.criticalEndRatios(i)*100;
                        fprintf(fid,'FRCfit(%s%1.1f%%) [l]\tLCIfit(%s%1.1f%%) [-]\t#BreathsFit(%s%1.1f%%) [-]\tdurationFit(%s%1.1f%%) [s]\tminuteVentilationFit(%s%1.1f%%) [ml/min]\tbreathRateFit(%s%1.1f%%) [1/min]\tdelta CEV Starndard to Fit(%s%1.1f%%) [l]\t', ...
                            gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio);
                        if momentRatio
                            fprintf(fid, 'mr0(%s%1.1f%%)\tmr1(%s%1.1f%%)\tmr2(%s%1.1f%%)\t', ...
                                         gasName,endRatio,gasName,endRatio,gasName,endRatio);
                        end
                    end
                end
                
                if TOevaluation
                    for i = 1:length(table.criticalTurnovers)
                        turnover=table.criticalTurnovers(i);
                        fprintf(fid,'FRCfromTO(%s%g) [l]\tCet(norm)fromTO(%s%g) [%%]\t#BreathsForTO(%s%g) [-]\tdurationForTO(%s%g) [s]\t', ...
                            gasName,turnover,gasName,turnover,gasName,turnover,gasName,turnover);
                    end
                end
            end
            fprintf(fid,'VolInspMean [ml]\tVolExpMean [ml]\trelVolInspStd [%%]\trelVolExpStd [%%]\t VolTidalMean [ml]\t');
        end
        
        if (n2mbw||tidal||dtgmbw) && capnoIndices
            fprintf(fid,'MeanFowlerDS [ml]\tMeanSnII_CO2 [1/L]\tMeanSnIII_CO2 [1/L]\tMeanfCO2E [%%]\tMeanCO2et [%%]\t');
            fprintf(fid,'StdFowlerDS [ml]\tStdSnII_CO2 [1/L]\tStdSnIII_CO2 [1/L]\tStdfCO2E [%%]\tStdCO2et [%%]\t');
            for i = 1:nBreathMax
                fprintf(fid,'FowlerDS(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'SnII_CO2(%d) [1/L]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'SnIII_CO2(%d) [1/L]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'fCO2E(%d) [%%]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'CO2et(%d) [%%]\t',i);
            end
            for i = 1:nBreathMax
                fprintf(fid,'includedCapnoBreaths(%d) [-]\t',i);
            end
        end
        
        if (n2mbw||tidal||dtgmbw) && tidalMean
            if table.TidalMeans.breath0 > 0
                modifier = sprintf('%d-%d',table.TidalMeans.startBreath,table.TidalMeans.endBreath);
            else
                modifier = sprintf('(%d)-(%d)',table.TidalMeans.breath0-table.TidalMeans.nBreath+1,table.TidalMeans.breath0);
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
        
        if (dtgsbw||n2mbw||tidal) && minVentState
            for i = 1:nBreathMax  %table.nHalfBreaths/2
                fprintf(fid,'inspiredVolume(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax  %table.nHalfBreaths/2
                fprintf(fid,'expiredVolume(%d) [ml]\t',i);
            end
            for i = 1:nBreathMax  %table.nHalfBreaths/2
                fprintf(fid,'meanInspiredFlow(%d) [ml/min]\t',i);
            end
            for i = 1:nBreathMax  %table.nHalfBreaths/2
                fprintf(fid,'meanExpiredFlow(%d) [ml/min]\t',i);
            end
        end
        
        if (dtgsbw||n2mbw||tidal||dtgmbw) && slopesState
            if torontoFile
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeCO2(%d) [%%/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeHe(%d) [%%/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeSF6(%d) [%%/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeCO2Norm(%d) [1/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeHeNorm(%d) [1/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeSF6Norm(%d) [1/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeHeNorm*VTexp(%d) [-]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeSF6Norm*VTexp(%d) [-]\t',i);
                end
            else
                if dtgsbw
                    for i = 1:nBreathMax  %table.nHalfBreaths/2
                        fprintf(fid,'SlopeCO2(%d) [%%/l]\t',i);
                    end
                    for i = 1:nBreathMax
                        fprintf(fid,'SlopeN2(%d) [%%/l]\t',i);
                    end
                else
                    for i = 1:nBreathMax  %table.nHalfBreaths/2
                        fprintf(fid,'SlopeCO2(%d) [%%/l]\t',i);
                    end
                    for i = 1:nBreathMax
                        fprintf(fid,'SlopeN2Norm(%d) [1/l]\t',i);
                    end
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeMMss(%d) [g/mol/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeMMssCalc(%d) [g/mol/l]\t',i);
                end
                for i = 1:nBreathMax
                    fprintf(fid,'SlopeDiffMMss(%d) [g/mol/l]\t',i);
                end
            end
        end
        
        if fastSlow                  
            for i = 1:length(table.N2.fastSlow)
                fields = fieldnames(table.N2.fastSlow{i});
                for j = 1:length(fields)
                    fprintf(fid, [fields{j} '\t']);
                end 
            end
            
            axis = {'N2 Concentration [%%]', 'Expired N2 Volume [L]'};
            for i = 1:length(axis)
                axisString = axis{i};
                for j = 1:60
                    fprintf(fid, [axisString ' Breath (' num2str(j) ')\t']);
                end 
            end
            %fprintf(fid, 'vTurbulent2\tvDiffusion2\tDiffusionFactor2 D\tR2\t');
            %fprintf(fid, 'vTurbulent3\tvLaminar3\tvDiffusion3\tDiffusionFactor3 D\tR2\t');
        end
        
        fprintf(fid,'\n');
    end

    calibration = parameters.Calibration;
    operation   = parameters.Operation;
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
    if torontoFile 
        fprintf(fid,'%s\t%s\t%s\t%s\t',datestr(now),nameInput,nameOutput,nameSettings);
    else
        fprintf(fid,'%s\t%s\t%s\t',datestr(now),nameInput,nameOutput);
        if n2mbw
            fprintf(fid,'%s\t',nameCalibration);
        end
        fprintf(fid,'%s\t',nameSettings);
        if n2mbw
            %fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t',calibration.delayFlow,calibration.coefficientsOptimal,operation.recalibration);
            fprintf(fid,'%e\t%d\t',calibration.delayFlow,operation.recalibration);
        end
    end
    fprintf(fid,'%d\t%d\t',operation.BTPS,simulation.staticState,simulation.dynamicState);
    fprintf(fid,'%e\t%e\t',calibration.factorBTPSinsp,calibration.factorBTPSexp);
    
    if dtgmbw
        gasByName=speciesNames{1};
        fprintf(fid, '%s\t', gasByName);
    end
    
    if dtgsbw
        fprintf(fid,'%e\t%e\t%e\t%e\t',table.CO2.slopePlain(1,2)*100/1000,table.N2.slopePlain(1,2)*100/1000,table.MMss.slope(1,2)*1000/1000,table.MMss.calcSlope(1,2)*1000/1000,table.MMss.diffSlope(1,2)*1000/1000);
        fprintf(fid,'%e\t%e\t',table.InterceptionVolume(1)*100,table.EndTidalDiffMMss(1)*1000,table.DiffMMssMax(1)*1000,table.DiffMMssMaxVolume(1)*100);
        fprintf(fid,'%e\t',table.volExp(1)*1000);
        fprintf(fid,'%e\t',table.minuteVentilationTotal(1)*1e6*60);
        fprintf(fid,'%e\t',table.breathRateTotal(1)*60);
    end
    
    if (n2mbw || dtgmbw) %&& ~torontoFile
        for j = 1 : nGases
            gas=table.(speciesNames{j});
            if ScondSacin
                data=gas.General.scondSacin;
                fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',data.Scond{1},data.Sacin{1},data.R2{1},length(data.InclusionIndices{1})+length(data.ExclusionIndices{1}),length(data.ExclusionIndices{1}));
                fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',data.Scond{2},data.Sacin{2},data.R2{2},length(data.InclusionIndices{2})+length(data.ExclusionIndices{2}),length(data.ExclusionIndices{2}));
            end
            if lciStandard
                data=gas.Standard;
                for i = 1:length(table.criticalEndRatios)
                    fprintf(fid,'%e\t%e\t%g\t%e\t%e\t%e\t',data.frc(i)*1e3,data.lci(i),data.nBreaths(i),data.duration(i),data.minuteVentilation(i)*1e6*60,data.breathRate(i)*60);
                    if momentRatio
                        fprintf(fid,'%e\t%e\t%g\t', gas.Standard.momentRatio(i,1), ...
                                                    gas.Standard.momentRatio(i,2), ...
                                                    gas.Standard.momentRatio(i,3));
                    end
                end
            end
            
            if lciByFit
                data=gas.ByFit;
                for i = 1:length(table.criticalEndRatios)
                    fprintf(fid,'%e\t%e\t%g\t%e\t%e\t%e\t%e\t', data.frc(i)*1e3,                            ...
                                                                data.lci(i),                                ...
                                                                data.nBreaths(i)-(1-data.partialBreath(i)), ...
                                                                data.duration(i),                           ...
                                                                data.minuteVentilation(i)*1e6*60,           ...
                                                                data.breathRate(i)*60,                      ...
                                                                data.deltaVolume(i)*1e3);
                    if momentRatio
                        fprintf(fid,'%e\t%e\t%g\t', gas.ByFit.momentRatio(i,1), ...
                                                    gas.ByFit.momentRatio(i,2), ...
                                                    gas.ByFit.momentRatio(i,3));
                    end
                end
            end
            
            if TOevaluation
                data=gas.Turnover;
                for i = 1:length(table.criticalTurnovers)
                    fprintf(fid,'%e\t%e\t%g\t%e\t',data.frc(i)*1e3,data.lci(i),data.nBreaths(i),data.duration(i));
                end
            end
        end
        
        data=table.Tidal;
        fprintf(fid,'%e\t%e\t%e\t',data.volInspMean*1e6,data.volExpMean*1e6,data.volInspStd*100,data.volExpStd*1e2, mean(table.volTidal)*1e6);
    end
    
    nHalfBreaths=floor((table.nHalfBreaths+1)/2);
    
    if (n2mbw||tidal|| dtgmbw) && capnoIndices
        capno=table.Capno;
        included=zeros(nHalfBreaths,1);
        included(capno.validIndices)=1;
        
        fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',capno.meanFowlerDead*1e6,capno.meanCoeffII/1e3,capno.meanCoeffIII/1e3,capno.meanfCO2E*1e2,capno.meanCO2et*1e2);
        fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',capno.stdFowlerDead*1e6, capno.stdCoeffII/1e3, capno.stdCoeffIII/1e3, capno.stdfCO2E*1e2, capno.stdCO2et*1e2 );
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.CO2.capno{j}.fowlerDead*1e6);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.CO2.capno{j}.coeffII/1e3);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.CO2.capno{j}.coeffIII/1e3);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.CO2.capno{j}.fCO2E*1e2);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.CO2.capno{j}.CO2et*1e2);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%g\t',included(j));
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
    end
    
    if (tidal || n2mbw || dtgmbw ) && tidalMean
        ref=table.TidalMeans;
        fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t',ref.minuteVentilation*1e6*60,ref.breathRate*60,ref.inspDuration,ref.expDuration,ref.timeInspPeakFlow,ref.timeExpPeakFlow,ref.volTidal*1e6,ref.tPTEF2tERatio);
    end
    
    if minVentState
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.volInsp(j)*1e6);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.volExp(j)*1e6);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.flowInspMean(j)*1e6*60);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
        for j = 1:min(nHalfBreaths,nBreathMax)
            fprintf(fid,'%e\t',table.flowExpMean(j)*1e6*60);
        end
        for j = nHalfBreaths+1:nBreathMax
            fprintf(fid,'\t');
        end
    end
    
    if dtgsbw
        scaleSpecies=100/1000;  % plain slopes
    else
        scaleSpecies=1/1000;    % slopes of normed species signals
    end
    
    if torontoFile
        if (n2mbw) && slopesState
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.CO2.slopePlain(i,2)*100/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.He.slopePlain(i,2)*100/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.SF6.slopePlain(i,2)*100/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.CO2.slope(i,2)*1/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.He.slope(i,2)*1/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.SF6.slope(i,2)*1/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.He.slope(i,2)*table.volExp(i));
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.SF6.slope(i,2)*table.volExp(i));
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end            
        end
    else
        if (dtgsbw||n2mbw||tidal||dtgmbw) && slopesState
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.CO2.slope(i,2)*100/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.N2.slope(i,2)*scaleSpecies);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.MMss.slope(i,2)*1000/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.MMss.calcSlope(i,2)*1000/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
            
            for i = 1:min(nHalfBreaths,nBreathMax)
                fprintf(fid,'%e\t',table.MMss.diffSlope(i,2)*1000/1000);
            end
            for i = nHalfBreaths+1:nBreathMax
                fprintf(fid,'\t');
            end
        end
    end
    
    if fastSlow 
        
        for i = 1:length(table.N2.fastSlow)
            
            fields = fieldnames(table.N2.fastSlow{i});
            for j = 1:length(fields)
                outputData = table.N2.fastSlow{i}.(fields{j});
                if isnumeric(outputData)
                   fprintf(fid,  '%d\t', outputData); 
                elseif ischar(outputData)
                    fprintf(fid, [outputData '\t']);
                else  
                   fprintf(fid, [func2str(outputData) '\t']);
                end
            end 
        end
        
        yCet = table.N2.General.cet_norm;                      
        yVN2 = (table.N2.exp-table.N2.reInsp)*1000;   
        axis = {yCet, yVN2};
        axisValues = zeros(60,1);
        
        for i = 1:length(axis)
            axisValues(1:length(axis{i}))     = axis{i}; 
            axisValues(length(axis{i})+1:end) = nan;
            for j = 1:60
                fprintf(fid, '%d\t', axisValues(j));
            end 
        end

%         fprintf(fid, '%d\t%d\t%d\t%d\t', ...
%                                          table.N2.twoCompartments.vTurbulent, ... 
%                                          table.N2.twoCompartments.vDiffusion, ...
%                                          table.N2.twoCompartments.D,          ...
%                                          table.N2.twoCompartments.r2);
%         fprintf(fid, '%d\t%d\t%d\t%d\t%d\t', ...
%                                          table.N2.threeCompartments.vTurbulent, ...
%                                          table.N2.threeCompartments.vLaminar,   ...
%                                          table.N2.threeCompartments.vDiffusion, ...
%                                          table.N2.threeCompartments.D,          ...
%                                          table.N2.threeCompartments.r2);
    end    

    fprintf(fid,'\n');
    fclose(fid);
    end
    if parameters.Simulation.reducedOutput
         if dtgsbw
        nBreathMax = 1; % to reduced the output column count (otherwise empty cells)
    end
    
    logFileName = 'logFile.txt';
    
    fid = fopen(logFileName, 'a');
    if fid<1
        error('cannot open file %s\n',logFileName);
    end
    
    finder=0;
    
    if ftell(fid)==0
        fprintf(fid,'Signal processing version\t');
        if torontoFile
%             fprintf(fid,'data+time\tinput file name\toutput file name\tsettings file name\t'); %detailed option
            fprintf(fid,'data+time\tinput file name\t'); %fast option
        else
%             fprintf(fid,'data+time\tinput file name\toutput file name\t'); %detailed option
            fprintf(fid,'data+time\tinput file name\t'); %fast option
            if n2mbw
%                 fprintf(fid,'calibration file name\t');
            end
%             fprintf(fid,'settings file name\t');
            if n2mbw
%                 fprintf(fid,'delay Iv->CO2\trecalibration\t');
                %fprintf(fid,'delay Iv->CO2\tdelay CO2->O2\tdelay CO2->MMss\ttauCO2\tkCO2\tkO2\tkMMss\trecalibration\t');
            end
        end
%         fprintf(fid,'BTPS\tstatic\tdynamic\t');
%         fprintf(fid,'BTPSinsp\tBTPSexp\t');
        if dtgmbw
            fprintf(fid,'Analysed Gas\t');
        end
        if dtgsbw
            fprintf(fid,'SlopeCO2 [%%/l]\tSlopeN2 [%%/l]\tSlopeMMss [g/mol/l]\tSlopeMMssCalc [g/mol/l]\tSlopeDiffMMss [g/mol/l]\t');
            fprintf(fid,'InterceptionVolume [%%/ExpVolume]\tEndTidalDiffMMss [g/mol]\tDiffMMssMax [g/mol]\tDiffMMssMaxVolume [%%ExpVol]\t');
            fprintf(fid,'ExpiredVolume [l]\t');
            fprintf(fid,'minuteVentilationManeuvre [ml/min]\t');
            fprintf(fid,'breathRateManeuvre [#/min]\t');
        end
        
        if n2mbw || dtgmbw
            for j = 1 : nGases
                if torontoFile
                    gasName=strcat('(',speciesNames{j},')');
                else
                    gasName='';
                end
                
                if ScondSacinNorm
                    normType = '';
                else
                    normType = '*VTexp';
                end
                
                if ScondSacin
                    if ScondSacinNorm
%                         fprintf(fid,'Scond%s [1/l/TO]\tSacin%s
%                         [1/l]\tR2%s [-]\t#Breaths total%s [-]\t#Breaths excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName); 
%                         fprintf(fid,'ScondStar%s [1/l/TO]\tSacinStar%s [1/l]\tR2%s [-]\t#BreathsStar total%s [-]\t#BreathsStar excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName);
                        fprintf(fid,'Scond%s [1/l/TO]\tSacin%s [1/l]\t',gasName,gasName);
                        fprintf(fid,'ScondStar%s [1/l/TO]\tSacinStar%s [1/l]\t' ,gasName,gasName);
                    else
%                         fprintf(fid,'Scond%s [1/TO]\tSacin%s [-]\tR2%s [-]\t#Breaths total%s [-]\t#Breaths excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName);
%                         fprintf(fid,'ScondStar%s [1/TO]\tSacinStar%s [-]\tR2%s [-]\t#BreathsStar total%s [-]\t#BreathsStar excluded%s [-]\t',gasName,gasName,gasName,gasName,gasName);
                        fprintf(fid,'Scond%s [1/l/TO]\tSacin%s [1/l]\t',gasName,gasName);
                        fprintf(fid,'ScondStar%s [1/l/TO]\tSacinStar%s [1/l]\t' ,gasName,gasName);
                    end
                end
                
                if torontoFile
                    gasName=strcat(speciesNames{j},':');
                else
                    gasName='';
                end
                if lciStandard
                    for i = 1:length(table.criticalEndRatios)
                        endRatio=table.criticalEndRatios(i)*100;
%                         fprintf(fid,'FRC(%s%1.1f%%) [l]\tLCI(%s%1.1f%%) [-]\t#Breaths(%s%1.1f%%) [-]\tduration(%s%1.1f%%) [s]\tminuteVentilation(%s%1.1f%%) [ml/min]\tbreathRate(%s%1.1f%%) [1/min]\t', ...
%                             gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio);
                         fprintf(fid,'FRC(%s%1.1f%%) [l]\tLCI(%s%1.1f%%) [-]\t', ...
                            gasName,endRatio,gasName,endRatio);
%                         if momentRatio
%                             fprintf(fid, 'mr0(%s%1.1f%%)\tmr1(%s%1.1f%%)\tmr2(%s%1.1f%%)\t', ...
%                                          gasName,endRatio,gasName,endRatio,gasName,endRatio);
%                         end
                    end
                end
                if lciByFit
                    for i = 1:length(table.criticalEndRatios)
                        endRatio=table.criticalEndRatios(i)*100;
%                         fprintf(fid,'FRCfit(%s%1.1f%%) [l]\tLCIfit(%s%1.1f%%) [-]\t#BreathsFit(%s%1.1f%%) [-]\tdurationFit(%s%1.1f%%) [s]\tminuteVentilationFit(%s%1.1f%%) [ml/min]\tbreathRateFit(%s%1.1f%%) [1/min]\tdelta CEV Starndard to Fit(%s%1.1f%%) [l]\t', ...
%                             gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio,gasName,endRatio);
                        fprintf(fid,'FRC(%s%1.1f%%) [l]\tLCI(%s%1.1f%%) [-]\t', ...
                            gasName,endRatio,gasName,endRatio);
%                         if momentRatio
%                             fprintf(fid, 'mr0(%s%1.1f%%)\tmr1(%s%1.1f%%)\tmr2(%s%1.1f%%)\t', ...
%                                          gasName,endRatio,gasName,endRatio,gasName,endRatio);
%                         end
                    end
                end
                
%                 if TOevaluation
%                     for i = 1:length(table.criticalTurnovers)
%                         turnover=table.criticalTurnovers(i);
%                         fprintf(fid,'FRCfromTO(%s%g) [l]\tCet(norm)fromTO(%s%g) [%%]\t#BreathsForTO(%s%g) [-]\tdurationForTO(%s%g) [s]\t', ...
%                             gasName,turnover,gasName,turnover,gasName,turnover,gasName,turnover);
%                     end
%                 end
            end
%             fprintf(fid,'VolInspMean [ml]\tVolExpMean [ml]\trelVolInspStd [%%]\trelVolExpStd [%%]\t VolTidalMean [ml]\t'); %detailed
            fprintf(fid,'VolTidalMean [ml]\t'); %fast
        end
        
%         if (n2mbw||tidal||dtgmbw) && capnoIndices
%             fprintf(fid,'MeanFowlerDS [ml]\tMeanSnII_CO2 [1/L]\tMeanSnIII_CO2 [1/L]\tMeanfCO2E [%%]\tMeanCO2et [%%]\t');
%             fprintf(fid,'StdFowlerDS [ml]\tStdSnII_CO2 [1/L]\tStdSnIII_CO2 [1/L]\tStdfCO2E [%%]\tStdCO2et [%%]\t');
%             for i = 1:nBreathMax
%                 fprintf(fid,'FowlerDS(%d) [ml]\t',i);
%             end
%             for i = 1:nBreathMax
%                 fprintf(fid,'SnII_CO2(%d) [1/L]\t',i);
%             end
%             for i = 1:nBreathMax
%                 fprintf(fid,'SnIII_CO2(%d) [1/L]\t',i);
%             end
%             for i = 1:nBreathMax
%                 fprintf(fid,'fCO2E(%d) [%%]\t',i);
%             end
%             for i = 1:nBreathMax
%                 fprintf(fid,'CO2et(%d) [%%]\t',i);
%             end
%             for i = 1:nBreathMax
%                 fprintf(fid,'includedCapnoBreaths(%d) [-]\t',i);
%             end
%         end
        
%         if (n2mbw||tidal||dtgmbw) && tidalMean
%             if table.TidalMeans.breath0 > 0
%                 modifier = sprintf('%d-%d',table.TidalMeans.startBreath,table.TidalMeans.endBreath);
%             else
%                 modifier = sprintf('(%d)-(%d)',table.TidalMeans.breath0-table.TidalMeans.nBreath+1,table.TidalMeans.breath0);
%             end
%             fprintf(fid,'MinuteVentilation%s [ml/min]\t',modifier);
%             fprintf(fid,'BreathRate%s [1/min]\t',modifier);
%             fprintf(fid,'inspDuration%s [s]\t',modifier);
%             fprintf(fid,'expDuration%s [s]\t',modifier);
%             fprintf(fid,'timeInspPeakFlow%s [s]\t',modifier);
%             fprintf(fid,'timeExpPeakFlow%s [s]\t',modifier);
%             fprintf(fid,'VolTidal%s [ml]\t',modifier);
%             fprintf(fid,'tPTEF2tERatio%s [-]\t',modifier);
%         end
        
%         if (dtgsbw||n2mbw||tidal) && minVentState
%             for i = 1:nBreathMax  %table.nHalfBreaths/2
%                 fprintf(fid,'inspiredVolume(%d) [ml]\t',i);
%             end
%             for i = 1:nBreathMax  %table.nHalfBreaths/2
%                 fprintf(fid,'expiredVolume(%d) [ml]\t',i);
%             end
%             for i = 1:nBreathMax  %table.nHalfBreaths/2
%                 fprintf(fid,'meanInspiredFlow(%d) [ml/min]\t',i);
%             end
%             for i = 1:nBreathMax  %table.nHalfBreaths/2
%                 fprintf(fid,'meanExpiredFlow(%d) [ml/min]\t',i);
%             end
%         end
        
%         if (dtgsbw||n2mbw||tidal||dtgmbw) && slopesState
%             if torontoFile
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeCO2(%d) [%%/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeHe(%d) [%%/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeSF6(%d) [%%/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeCO2Norm(%d) [1/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeHeNorm(%d) [1/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeSF6Norm(%d) [1/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeHeNorm*VTexp(%d) [-]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeSF6Norm*VTexp(%d) [-]\t',i);
%                 end
%             else
%                 if dtgsbw
%                     for i = 1:nBreathMax  %table.nHalfBreaths/2
%                         fprintf(fid,'SlopeCO2(%d) [%%/l]\t',i);
%                     end
%                     for i = 1:nBreathMax
%                         fprintf(fid,'SlopeN2(%d) [%%/l]\t',i);
%                     end
%                 else
%                     for i = 1:nBreathMax  %table.nHalfBreaths/2
%                         fprintf(fid,'SlopeCO2(%d) [%%/l]\t',i);
%                     end
%                     for i = 1:nBreathMax
%                         fprintf(fid,'SlopeN2Norm(%d) [1/l]\t',i);
%                     end
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeMMss(%d) [g/mol/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeMMssCalc(%d) [g/mol/l]\t',i);
%                 end
%                 for i = 1:nBreathMax
%                     fprintf(fid,'SlopeDiffMMss(%d) [g/mol/l]\t',i);
%                 end
%             end
%         end
        
%         if fastSlow                  
%             for i = 1:length(table.N2.fastSlow)
%                 fields = fieldnames(table.N2.fastSlow{i});
%                 for j = 1:length(fields)
%                     fprintf(fid, [fields{j} '\t']);
%                 end 
%             end
%             
%             axis = {'N2 Concentration [%%]', 'Expired N2 Volume [L]'};
%             for i = 1:length(axis)
%                 axisString = axis{i};
%                 for j = 1:60
%                     fprintf(fid, [axisString ' Breath (' num2str(j) ')\t']);
%                 end 
%             end
%             %fprintf(fid, 'vTurbulent2\tvDiffusion2\tDiffusionFactor2 D\tR2\t');
%             %fprintf(fid, 'vTurbulent3\tvLaminar3\tvDiffusion3\tDiffusionFactor3 D\tR2\t');
%         end
        
        fprintf(fid,'\n');
    end
%________________________________________________________________________________________________________________________________________
finder=0;

    calibration = parameters.Calibration;
    operation   = parameters.Operation;
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
    if torontoFile 
%         fprintf(fid,'%s\t%s\t%s\t%s\t',datestr(now),nameInput,nameOutput,nameSettings); %detailed
        fprintf(fid,'%s\t%s\t',datestr(now),nameInput); %fast
    else
%         fprintf(fid,'%s\t%s\t%s\t',datestr(now),nameInput,nameOutput);
        fprintf(fid,'%s\t%s\t',datestr(now),nameInput);
        if n2mbw
%             fprintf(fid,'%s\t',nameCalibration);
        end
%         fprintf(fid,'%s\t',nameSettings);
        if n2mbw
            %fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%d\t',calibration.delayFlow,calibration.coefficientsOptimal,operation.recalibration);
%             fprintf(fid,'%e\t%d\t',calibration.delayFlow,operation.recalibration);
        end
    end
%     fprintf(fid,'%d\t%d\t',operation.BTPS,simulation.staticState,simulation.dynamicState);
%     fprintf(fid,'%e\t%e\t',calibration.factorBTPSinsp,calibration.factorBTPSexp);
    
    if dtgmbw
        gasByName=speciesNames{1};
        fprintf(fid, '%s\t', gasByName);
    end
    
    if dtgsbw
        fprintf(fid,'%e\t%e\t%e\t%e\t',table.CO2.slopePlain(1,2)*100/1000,table.N2.slopePlain(1,2)*100/1000,table.MMss.slope(1,2)*1000/1000,table.MMss.calcSlope(1,2)*1000/1000,table.MMss.diffSlope(1,2)*1000/1000);
        fprintf(fid,'%e\t%e\t',table.InterceptionVolume(1)*100,table.EndTidalDiffMMss(1)*1000,table.DiffMMssMax(1)*1000,table.DiffMMssMaxVolume(1)*100);
        fprintf(fid,'%e\t',table.volExp(1)*1000);
        fprintf(fid,'%e\t',table.minuteVentilationTotal(1)*1e6*60);
        fprintf(fid,'%e\t',table.breathRateTotal(1)*60);
    end
    
    if (n2mbw || dtgmbw) %&& ~torontoFile
        for j = 1 : nGases
            gas=table.(speciesNames{j});
            if ScondSacin
                data=gas.General.scondSacin;
%                 fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',data.Scond{1},data.Sacin{1},data.R2{1},length(data.InclusionIndices{1})+length(data.ExclusionIndices{1}),length(data.ExclusionIndices{1}));
%                 fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',data.Scond{2},data.Sacin{2},data.R2{2},length(data.InclusionIndices{2})+length(data.ExclusionIndices{2}),length(data.ExclusionIndices{2}));
                fprintf(fid,'%g\t%g\t',data.Scond{1},data.Sacin{1});
                fprintf(fid,'%g\t%g\t',data.Scond{2},data.Sacin{2});

            
            end
            if lciStandard
                data=gas.Standard;
                for i = 1:length(table.criticalEndRatios)
%                     fprintf(fid,'%e\t%e\t%g\t%e\t%e\t%e\t',data.frc(i)*1e3,data.lci(i),data.nBreaths(i),data.duration(i),data.minuteVentilation(i)*1e6*60,data.breathRate(i)*60);
                    fprintf(fid,'%e\t%e\t',data.frc(i)*1e3,data.lci(i));
%                     if momentRatio
%                         fprintf(fid,'%e\t%e\t%g\t', gas.Standard.momentRatio(i,1), ...
%                                                     gas.Standard.momentRatio(i,2), ...
%                                                     gas.Standard.momentRatio(i,3));
%                     end
                end
            end
            
            if lciByFit
                data=gas.ByFit;
                for i = 1:length(table.criticalEndRatios)
%                     fprintf(fid,'%e\t%e\t%g\t%e\t%e\t%e\t%e\t', data.frc(i)*1e3,                            ...
%                                                                 data.lci(i),                                ...
%                                                                 data.nBreaths(i)-(1-data.partialBreath(i)), ...
%                                                                 data.duration(i),                           ...
%                                                                 data.minuteVentilation(i)*1e6*60,           ...
%                                                                 data.breathRate(i)*60,                      ...
%                                                                 data.deltaVolume(i)*1e3);
                     fprintf(fid,'%e\t%e\t', data.frc(i)*1e3, data.lci(i));
%                     if momentRatio
%                         fprintf(fid,'%e\t%e\t%g\t', gas.ByFit.momentRatio(i,1), ...
%                                                     gas.ByFit.momentRatio(i,2), ...
%                                                     gas.ByFit.momentRatio(i,3));
%                     end
                end
            end
            
%             if TOevaluation
%                 data=gas.Turnover;
%                 for i = 1:length(table.criticalTurnovers)
%                     fprintf(fid,'%e\t%e\t%g\t%e\t',data.frc(i)*1e3,data.lci(i),data.nBreaths(i),data.duration(i));
%                 end
%             end
        end
        
        data=table.Tidal;
        fprintf(fid,'%e\t', mean(table.volTidal)*1e6);
    end
    
%     nHalfBreaths=floor((table.nHalfBreaths+1)/2);
    
%     if (n2mbw||tidal|| dtgmbw) && capnoIndices
%         capno=table.Capno;
%         included=zeros(nHalfBreaths,1);
%         included(capno.validIndices)=1;
%         
%         fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',capno.meanFowlerDead*1e6,capno.meanCoeffII/1e3,capno.meanCoeffIII/1e3,capno.meanfCO2E*1e2,capno.meanCO2et*1e2);
%         fprintf(fid,'%g\t%g\t%g\t%g\t%g\t',capno.stdFowlerDead*1e6, capno.stdCoeffII/1e3, capno.stdCoeffIII/1e3, capno.stdfCO2E*1e2, capno.stdCO2et*1e2 );
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.CO2.capno{j}.fowlerDead*1e6);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.CO2.capno{j}.coeffII/1e3);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.CO2.capno{j}.coeffIII/1e3);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.CO2.capno{j}.fCO2E*1e2);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.CO2.capno{j}.CO2et*1e2);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%g\t',included(j));
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%     end
    
%     if (tidal || n2mbw || dtgmbw ) && tidalMean
%         ref=table.TidalMeans;
%         fprintf(fid,'%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t',ref.minuteVentilation*1e6*60,ref.breathRate*60,ref.inspDuration,ref.expDuration,ref.timeInspPeakFlow,ref.timeExpPeakFlow,ref.volTidal*1e6,ref.tPTEF2tERatio);
%     end
    
%     if minVentState
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.volInsp(j)*1e6);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.volExp(j)*1e6);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.flowInspMean(j)*1e6*60);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%         for j = 1:min(nHalfBreaths,nBreathMax)
%             fprintf(fid,'%e\t',table.flowExpMean(j)*1e6*60);
%         end
%         for j = nHalfBreaths+1:nBreathMax
%             fprintf(fid,'\t');
%         end
%     end
    
%     if dtgsbw
%         scaleSpecies=100/1000;  % plain slopes
%     else
%         scaleSpecies=1/1000;    % slopes of normed species signals
%     end
    
%     if torontoFile
%         if (n2mbw) && slopesState
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.CO2.slopePlain(i,2)*100/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.He.slopePlain(i,2)*100/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.SF6.slopePlain(i,2)*100/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.CO2.slope(i,2)*1/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.He.slope(i,2)*1/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.SF6.slope(i,2)*1/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.He.slope(i,2)*table.volExp(i));
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.SF6.slope(i,2)*table.volExp(i));
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end            
%         end
%     else
%         if (dtgsbw||n2mbw||tidal||dtgmbw) && slopesState
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.CO2.slope(i,2)*100/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.N2.slope(i,2)*scaleSpecies);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.MMss.slope(i,2)*1000/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.MMss.calcSlope(i,2)*1000/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%             
%             for i = 1:min(nHalfBreaths,nBreathMax)
%                 fprintf(fid,'%e\t',table.MMss.diffSlope(i,2)*1000/1000);
%             end
%             for i = nHalfBreaths+1:nBreathMax
%                 fprintf(fid,'\t');
%             end
%         end
%     end
    
%     if fastSlow 
%         
%         for i = 1:length(table.N2.fastSlow)
%             
%             fields = fieldnames(table.N2.fastSlow{i});
%             for j = 1:length(fields)
%                 outputData = table.N2.fastSlow{i}.(fields{j});
%                 if isnumeric(outputData)
%                    fprintf(fid,  '%d\t', outputData); 
%                 elseif ischar(outputData)
%                     fprintf(fid, [outputData '\t']);
%                 else  
%                    fprintf(fid, [func2str(outputData) '\t']);
%                 end
%             end 
%         end
%         
%         yCet = table.N2.General.cet_norm;                      
%         yVN2 = (table.N2.exp-table.N2.reInsp)*1000;   
%         axis = {yCet, yVN2};
%         axisValues = zeros(60,1);
%         
%         for i = 1:length(axis)
%             axisValues(1:length(axis{i}))     = axis{i}; 
%             axisValues(length(axis{i})+1:end) = nan;
%             for j = 1:60
%                 fprintf(fid, '%d\t', axisValues(j));
%             end 
%         end
% 
% %         fprintf(fid, '%d\t%d\t%d\t%d\t', ...
% %                                          table.N2.twoCompartments.vTurbulent, ... 
% %                                          table.N2.twoCompartments.vDiffusion, ...
% %                                          table.N2.twoCompartments.D,          ...
% %                                          table.N2.twoCompartments.r2);
% %         fprintf(fid, '%d\t%d\t%d\t%d\t%d\t', ...
% %                                          table.N2.threeCompartments.vTurbulent, ...
% %                                          table.N2.threeCompartments.vLaminar,   ...
% %                                          table.N2.threeCompartments.vDiffusion, ...
% %                                          table.N2.threeCompartments.D,          ...
% %                                          table.N2.threeCompartments.r2);
%     end    

    fprintf(fid,'\n');
    fclose(fid);
    end
    ok=1;
    
end

