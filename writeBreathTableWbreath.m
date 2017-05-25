%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to write ASCII breath table data
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.3, 30. Jan. 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok,table] = writeBreathTableWbreath(table,file,parameters,flag)

    SF6             =   parameters.Simulation.SF6;
    
    indexCrit       =   parameters.Operation.nVentilation;      % #breaths for minute ventilation evaluation
    TOevaluation    =   parameters.Simulation.TOevaluation;
    lciStandard    	=   parameters.Simulation.lciStandard;
    lciByFit      	=   parameters.Simulation.lciByFit;
    
    deltaState      =   parameters.Simulation.deltaMMState;
    
    tidalMean       =   parameters.Simulation.tidalMeanState;
    
    cropping        =   parameters.Simulation.croppingState;	% data cropping
    washout         =   parameters.Simulation.washout;          % washout state
    tidal           =   parameters.Simulation.tidalState;       % tidal state
    
    DSpost          =   parameters.Device.volumePostcap;
    DSpre           =   parameters.Device.volumePrecap;
    DS              =   DSpost+DSpre;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Auto-calibration of limit values for washin / washout
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    minMMInsp       =   min(table.MMInspMean);
    maxMMInsp       =   max(table.MMInspMean);
    
    hysteresisInsp  =   (maxMMInsp-minMMInsp)/3;
    
    if SF6
        MMcritWI        =   minMMInsp+0.001;
        MMcritWO        =   maxMMInsp-0.001;
    else
        MMcritWO        =   minMMInsp+hysteresisInsp;
        MMcritWI        =   maxMMInsp-hysteresisInsp;
    end
    
    tidalOnly       =   0;
    if maxMMInsp-minMMInsp < 0.002
        tidalOnly   =   1;
    end
    
    
    % old discrimination limit values
    % MMcritWI        =   parameters.Operation.MMair+0.002;
    % MMcritWO        =   parameters.Operation.MMTG -0.002;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation of tidal breath table
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if tidal
        if tidalOnly==1
            preBreaths=find(table.MMInspMean>0);
            wiBreaths=[];
            woBreaths=[];
        else
            wiBreaths=find(table.MMInspMean>MMcritWI);
            woBreaths=find(table.MMInspMean<MMcritWO);
            boundary=find(diff(woBreaths)>1,1);
            if isempty(boundary)
                if woBreaths(1)==1
                    preBreaths=woBreaths;
                    woBreaths=[];
                else
                    preBreaths=wiBreaths;
                    wiBreaths=[];
                end
            else
                preBreaths=woBreaths(1:boundary);
                woBreaths=woBreaths(boundary+1:end);
            end
        end
        table.breathType=cell(length(table.MMInspMean),1);
        table.breathType(preBreaths)={'pre'};
        table.breathType(woBreaths) ={'wo'};
        table.breathType(wiBreaths) ={'wi'};
        table.DeeBase=0;    % dummy value
        
        table.CEV=cumsum(table.VolExp);
        timeLine=cumsum(table.breathDuration);
        table.minuteVentilationTotal=table.CEV(end,1)/timeLine(end,1);
        table.breathRateTotal=length(table.CEV(:,1))/timeLine(end,1);
        if indexCrit<=length(table.CEV(:,1))
            table.minuteVentilationCrit=table.CEV(indexCrit,1)/timeLine(indexCrit,1);
            table.breathRateCrit=indexCrit/timeLine(indexCrit,1);
        else
            fprintf('WARNING: cannot determine minute ventilation, not enough breaths\n');
            table.minuteVentilationCrit=0;
            table.breathRateCrit=0;
        end
        
        [ok,table] = tidalAnalysis(table,parameters);
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation of MBW breath table, autodetection of washin or washout
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Detection of MBW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if SF6 == 1
        mbwBreaths=find(table.MMInspMean>MMcritWI);
    else
        mbwBreaths=find(table.MMInspMean<MMcritWO);
    end
    if isempty(mbwBreaths)
        error('writeBreathTableWbreath:noConsistentData','(writeBreathTable): cannot determine MBW breaths (MM !> %g g/mol) during inspiration)',MMcritWI*1000);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Overriding washin/washout state due to user selection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if cropping
        if length(mbwBreaths)>length(table.VolMMInsp)-mbwBreaths(end)
            washout=0;
            fprintf('manover is detected as a washin\n')
        else
            washout=1;
            fprintf('manover is detected as a washout\n')
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of washin/washout related MM parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if washout==1	% is washout
        if SF6 == 1
            mbwBreaths  =   find(table.MMInspMean<MMcritWO);
        else
            mbwBreaths  =   find(table.MMInspMean>MMcritWO);
        end
        if isempty(mbwBreaths)
            error('writeBreathTableWbreath:notEnoughPreBreaths','(writeBreathTable): cannot determine MBWO breaths (MM !< %g g/mol) during inspiration)',MMcritWO*1000);
        end
        lastPart        =   find(diff(mbwBreaths)>1,1,'last');
        if ~isempty(lastPart)
            mbwBreaths  =  	mbwBreaths(lastPart+1:end);
        end
        
        startIndex          =   mbwBreaths(1);
        endIndex            =   mbwBreaths(end);
        table.CeiEnd        =   mean(table.MMei(endIndex-2:endIndex,1));
        table.MMEnd         =   mean(table.MMExpMean(endIndex-2:endIndex,1));
        table.VolMMExpEff   =   table.VolMMExp-table.VolExp.*table.MMEnd;
        if SF6 == 1
            table.VolMMExpEff	=   table.VolMMExp-table.VolExp.*table.MMEnd;
        else
            table.VolMMExpEff	=   -(table.VolMMExp-table.VolExp.*table.MMEnd);
        end
    elseif washout==0	% is washin
        if SF6 == 1
            mbwBreaths      =   find(table.MMInspMean>MMcritWI);
        else
            mbwBreaths      =   find(table.MMInspMean<MMcritWI);
        end
        startIndex          =   mbwBreaths(1);
        endIndex            =   mbwBreaths(end);
        table.CeiEnd        =   mean(table.MMei(endIndex-2:endIndex,1));
        table.MMEnd         =   mean(table.MMExpMean(endIndex-2:endIndex,1));
        table.VolMMExpEff   =  -(table.VolMMExp-table.VolExp.*table.MMEnd);
        if SF6 == 1
            table.VolMMExpEff	=   -(table.VolMMExp-table.VolExp.*table.MMEnd);
        else
            table.VolMMExpEff	=   table.VolMMExp-table.VolExp.*table.MMEnd;
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Setting mbwBreaths and running wBreath algorithm
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mbwBreaths=(startIndex:endIndex);
    
    if startIndex<4
        error('writeBreathTableWbreath:notEnoughPreBreaths','(writeBreathTable): cannot determine CeeStart, not enougth pre-breaths (%d)',startIndex-1);
    else
        table.CeeStart      =   mean(table.MMee(startIndex-3:startIndex-1,1));
        table.CeiStart      =   mean(table.MMei(startIndex-3:startIndex-1,1));
        table.DeeBaseFRC	=   table.CeiEnd-table.CeiStart;
        table.CeeEndFRC     =   table.CeeStart+table.DeeBaseFRC;
        table.DeeFRC        =   table.CeeEndFRC-table.MMee;
        table.DeeNormFRC    =   table.DeeFRC/table.DeeBaseFRC;
        if deltaState
            table.DeeBaseLCI	=   sign(table.DeeBaseFRC)*parameters.Simulation.deltaMMFixed;
        else
            table.DeeBaseLCI	=   table.DeeBaseFRC;
        end
        table.CeeEndLCI     =   table.CeeStart+table.DeeBaseLCI;
        table.DeeLCI        =   table.CeeEndLCI-table.MMee;
        table.DeeNormLCI    =   table.DeeLCI/table.DeeBaseLCI;
        minDeeNormLCI       =   min(table.DeeNormLCI);
        if minDeeNormLCI < 0
            fprintf('WARNING (writeBreathTableWbreath): data inconsistent, negative tracer concentations!\n')
        end
        stabilityPreBreaths =   std(table.MMei(startIndex-3:startIndex-1,1))/abs(table.DeeBaseLCI);
        if stabilityPreBreaths>0.05
            error('writeBreathTableWbreath:dataInconsistent','(writeBreathTable): too large variation in pre breaths MMei (%g)',stabilityPreBreaths)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Cropping table to MBW breaths only
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fields=fieldnames(table);
    for i=1:numel(fields)
        if length(table.(fields{i}))>length(mbwBreaths)
            %         disp(fields{i});
            table.(fields{i})=table.(fields{i})(mbwBreaths,:);
        end
    end
    table.nHalfBreaths=2*length(mbwBreaths);
    
    if ~isempty(find(table.DeeLCI<0, 1))
        table.DeeFRC        =   abs(table.DeeFRC);
        table.DeeNormFRC    =   abs(table.DeeNormFRC);
        table.DeeLCI        =   abs(table.DeeLCI);
        table.DeeNormLCI    =   abs(table.DeeNormLCI);
        fprintf('WARNING (writeBreathTableWbreath): data inconsistent as negative concentrations (DeeLCI) found!\n')
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculating extra fields, FRC, LCI, etc.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    table.tracerVol         =   table.VolMMExpEff/abs(table.DeeBaseFRC);
    table.tracerVolReInsp   =   DSpost*table.DeeNormFRC;
    table.tracerVolNetto    =   table.tracerVol-table.tracerVolReInsp;
    table.cumTracerVolNetto =   cumsum(table.tracerVolNetto);
    
    table.CEV               =   cumsum(table.VolExp);
    table.CEV_DS            =   table.CEV-(1:length(mbwBreaths))'*DS;
    table.FRCSP             =   table.cumTracerVolNetto./(1-table.DeeNormFRC);
    table.FRCAO             =   table.FRCSP-DSpre;
    table.TO                =   table.CEV_DS./table.FRCAO;
    table.minuteVentilation =   table.VolExp./table.breathDuration;
    
    fprintf('DeeBaseFRC = %f g/mol\n', table.DeeBaseFRC*1000);
    fprintf('DeeBaseLCI = %f g/mol\n', table.DeeBaseLCI*1000);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % User output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % % % % % table.VolInspMean=mean(table.VolInsp);
    % % % % % table.VolExpMean =mean(table.VolExp );
    % % % % % table.VolInspStd =std(table.VolInsp,1)/table.VolInspMean;
    % % % % % table.VolExpStd  =std(table.VolExp ,1)/table.VolExpMean;
    fprintf('mean(VolInsp)=%.1f ml+-%.1f%%, mean(VolExp)=%.1f ml+-%.1f%%\n', table.VolInspMean*1e6,table.VolInspStd*100,table.VolExpMean*1e6,table.VolExpStd*100);
    
    table.criticalEndRatios=[0.01,0.5,1.0,1.5,2.0,2.5,3,4,5,7.5,10,15,20]/100;  	% critical end values
    table.criticalTurnovers=[6,8,10];                                               % critical turnover values
    
    if lciStandard
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detection of MBW end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        table.targetIndexConcentration=zeros(size(table.criticalEndRatios));        % declaration of targetIndices array
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Dee(norm) based determination of TO=LCI
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:length(table.criticalEndRatios)
            targetIndices=find(table.DeeNormLCI(:,1)<table.criticalEndRatios(i));
            lciIndex=getLCIIndex(targetIndices,0);
            if lciIndex>0
                table.targetIndexConcentration(i)=lciIndex;
            end
        end
        
        Direct.DeeNorm=zeros(size(table.criticalEndRatios));
        Direct.LCI=zeros(size(table.criticalEndRatios));
        Direct.FRC=zeros(size(table.criticalEndRatios));
        Direct.duration=zeros(size(table.criticalEndRatios));
        Direct.minuteVentilation=zeros(size(table.criticalEndRatios));
        Direct.breathRate=zeros(size(table.criticalEndRatios));
        Direct.nBreaths=zeros(size(table.criticalEndRatios));
        
        disp(' ');
        disp('FRC and LCI for given critical Dee');
        for i = 1:length(table.criticalEndRatios)
            DeeCrit=table.criticalEndRatios(i);
            if table.targetIndexConcentration(i)>0
                index=table.targetIndexConcentration(i);
                Direct.DeeNorm(i)=DeeCrit;
                Direct.LCI(i)=table.TO(index,1);
                Direct.FRC(i)=table.FRCAO(index,1);
                
                [m0,m1,m2]=moments(table.CEV_DS/Direct.FRC(i),table.DeeNormLCI,Direct.LCI(i));
                Direct.m0(i)=m0;
                Direct.m1(i)=m1;
                Direct.m2(i)=m2;
                Direct.r1(i)=m1/m0;
                Direct.r2(i)=m2/m0;
                
                Direct.duration(i)=sum(table.breathDuration(1:index,1));
                Direct.minuteVentilation(i)=table.CEV(index,1)/Direct.duration(i);
                Direct.nBreaths(i)=index;
                Direct.breathRate(i)=index/Direct.duration(i);
                fprintf('#breath=%3g: duration=%5.1f s, FRC=%6.1f ml, LCI=%5.2f, MR1=%5.2f, MR2=%5.2f at Dee(norm)=%4.1f %%, BR=%4.1f 1/min, MV=%6.1f ml/min\n', ...
                    Direct.nBreaths(i),Direct.duration(i),Direct.FRC(i)*1e6,Direct.LCI(i),Direct.r1(i),Direct.r2(i),DeeCrit*100,Direct.breathRate(i)*60,Direct.minuteVentilation(i)*1e6*60);
            else
                fprintf('WARNING: cannot determine FRC and LCI for Dee(norm) = %f %% (MM stays above thershold)\n',DeeCrit*100);
            end
        end
        table.Direct = Direct;
    end
    
    timeLine=cumsum(table.breathDuration);
    table.minuteVentilationTotal=table.CEV(end,1)/timeLine(end,1);
    table.breathRateTotal=length(table.CEV(:,1))/timeLine(end,1);
    if indexCrit<=length(table.CEV(:,1))
        table.minuteVentilationCrit=table.CEV(indexCrit,1)/timeLine(indexCrit,1);
        table.breathRateCrit=indexCrit/timeLine(indexCrit,1);
    else
        fprintf('WARNING: cannot determine minute ventilation, not enought breaths\n');
        table.minuteVentilationCrit=0;
        table.breathRateCrit=0;
    end
    
    if TOevaluation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Turnover based detection of Dee(norm)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Turnover.DeeNorm=zeros(size(table.criticalTurnovers));
        Turnover.LCI=zeros(size(table.criticalTurnovers));
        Turnover.FRC=zeros(size(table.criticalTurnovers));
        Turnover.duration=zeros(size(table.criticalTurnovers));
        Turnover.nBreaths=zeros(size(table.criticalTurnovers));
        
        disp(' ');
        disp('FRC and Dee(norm) for given TO');
        for i = 1:length(table.criticalTurnovers)
            TOcrit = table.criticalTurnovers(i);
            targetIndices=find(table.TO(:,1)>TOcrit);
            if ~isempty(targetIndices)
                index = targetIndices(1);    % due to cropping, the startindex is now generally equal to 1
                Turnover.DeeNorm(i)=table.DeeNormLCI(index);
                Turnover.LCI(i)=TOcrit;
                Turnover.FRC(i)=table.FRCAO(index);
                
                [m0,m1,m2]=moments(table.CEV_DS/Turnover.FRC(i),table.DeeNormLCI,Turnover.LCI(i));
                Turnover.m0(i)=m0;
                Turnover.m1(i)=m1;
                Turnover.m2(i)=m2;
                Turnover.r1(i)=m1/m0;
                Turnover.r2(i)=m2/m0;
                
                Turnover.duration(i)=sum(table.breathDuration(1:index,1));
                Turnover.nBreaths(i)=index;
                fprintf('#breath=%d: duration=%.1f s, FRC=%.1f ml, MR1=%5.2f, MR2=%5.2f, Dee(norm)=%.2f %% at TO=%g\n', ...
                    Turnover.nBreaths(i),Turnover.duration(i),Turnover.FRC(i)*1e6,Turnover.r1(i),Turnover.r2(i),Turnover.DeeNorm(i)*100,TOcrit);
            else
                Turnover.m0(i)=0;
                Turnover.m1(i)=0;
                Turnover.m2(i)=0;
                Turnover.r1(i)=0;
                Turnover.r2(i)=0;
                Turnover.duration(i)=0;
                Turnover.nBreaths(i)=0;
                fprintf('WARNING: cannot determine Dee(norm) and FRC for TO=%g (TOmax<TOcrit)\n',TOcrit);
            end
        end
        table.Turnover = Turnover;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fitting linear slopes in semilog plot (used only as a guide for the eye in graphs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    table.logDeeNorm = log(table.DeeNormLCI);
    CEVmax = max(table.CEV_DS);
    [fitParametersUpper,fitFunctionUpper,cevUpper,fitUpper,rmsUpper,r2Upper]=getFit('lin',0,         CEVmax/3,table.CEV_DS,table.logDeeNorm);
    [fitParametersLower,fitFunctionUpper,cevLower,fitLower,rmsLower,r2Lower]=getFit('lin',2*CEVmax/3,CEVmax,  table.CEV_DS,table.logDeeNorm);
    table.DeeNormFitUpper = exp(fitUpper);
    table.DeeNormFitLower = exp(fitLower);
    table.CEVUpper = cevUpper;
    table.CEVLower = cevLower;
    
    if lciByFit
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % LCI = LCI(Dee_norm) determination using double exponential fits
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Fit.DeeNorm=zeros(size(table.criticalEndRatios));
        Fit.LCI=zeros(size(table.criticalEndRatios));
        Fit.FRC=zeros(size(table.criticalEndRatios));
        Fit.duration=zeros(size(table.criticalEndRatios));
        Fit.nBreaths=zeros(size(table.criticalEndRatios));
        
        disp(' ');
        disp('LCI calculation by fitting Dee(norm) and FRC');
        for i = 1:length(table.criticalEndRatios)
            DeeCrit = table.criticalEndRatios(i);
            lastIndex=min(min(find(table.DeeNormLCI<DeeCrit,1,'first'))+2,length(table.DeeNormLCI));
            if ~isempty(lastIndex)
                CEVmax = table.CEV_DS(lastIndex);
                CEVmin = CEVmax*0.5;	% arbitrarily chosen, used as lower bound for linear fit of FRC
                [parametersDee,handleDee,cevDee,fitDee,rmsDee,r2Cee]=getFit('nonlin',0,CEVmax,table.CEV_DS,table.DeeNormLCI);
                [parametersFRC,handleFRC,cevFRC,fitFRC,rmsFRC,r2FRC]=getFit('lin',CEVmin,CEVmax,table.CEV_DS,table.FRCAO);
                [parametersDuration,handleDuration,cevDuration,fitDuration,rmsDuration,r2Duration]=getFit('lin',0,CEVmax,table.CEV_DS,cumsum(table.breathDuration));
                
                if handleDee(CEVmax*20)<DeeCrit
                    CEV = fzero(@(v) handleDee(v)-DeeCrit,CEVmax*2);  % finding cumulative expired volume at CeeCrit
                else
                    CEV = nan;
                end
            end
            if ~isempty(lastIndex) && ~isnan(CEV)
                FRC = handleFRC(CEV);
                LCI = CEV/FRC;
                duration=handleDuration(CEV);
                Fit.DeeNorm(i)=DeeCrit;
                Fit.LCI(i)=LCI;
                Fit.FRC(i)=FRC;
                
                [m0,m1,m2]=moments(table.CEV_DS/Fit.FRC(i),table.DeeNormLCI,Fit.LCI(i));
                Fit.m0(i)=m0;
                Fit.m1(i)=m1;
                Fit.m2(i)=m2;
                Fit.r1(i)=m1/m0;
                Fit.r2(i)=m2/m0;
                
                Fit.nBreaths(i)=find(table.DeeNormLCI<DeeCrit,1,'first');
                Fit.duration(i)=duration;
                Fit.handleDee{i}=handleDee;
                Fit.handleFRC{i}=handleFRC;
                Fit.minuteVentilation(i)=CEV/duration;
                Fit.breathRate(i)=Fit.nBreaths(i)/duration;
                fprintf('#breath=%3g: duration=%5.1f s, FRC=%6.1f ml, LCI=%5.2f, MR1=%5.2f, MR2=%5.2f at Dee(norm)=%4.1f %%, BR=%4.1f 1/min, MV=%6.1f ml/min\n', ...
                    Fit.nBreaths(i),Fit.duration(i),Fit.FRC(i)*1e6,Fit.LCI(i),Fit.r1(i),Fit.r2(i),DeeCrit*100,Fit.breathRate(i)*60,Fit.minuteVentilation(i)*1e6*60);
            else
                fprintf('WARNING: cannot determine FRC and LCI for Dee = %f %% (MM stays above thershold)\n',DeeCrit*100);
            end
        end
        table.Fit = Fit;
    end
    
    if tidalMean
        [ok,table] = tidalAnalysis(table,parameters);
    end
    
    
    fid = fopen(file, 'w');
    if fid<1
        error('writeBreathTableWbreath:cannotOpenFile','cannot open file %s',file)
    end
    
    header1='Breath #\tDeeFRC [%%]\tDeeLCI [%%]\tTO [-]\tFRC [l]\tDeeFRC(norm) [%%]\tDeeLCI(norm) [%%]\tVolInsp [l]\tVolExp [l]\t'; % 9
    header2='CEV [l]\tCEV-DS [l]\ttracerVolExp [ml]\ttracerVolNetto [ml]\tCumTracerVolNetto [ml]\ttracerVolReInsp [ml]\t'; % 6
    header3='FlowInsp.mean [ml/s]\tFlowExp.mean [ml/s]\t'; % 2
    header = strcat(header1,header2,header3,'\r\n');
    format='';
    nColumns=17; % = 9+6+2
    for i=1:nColumns
        format=strcat(format,'%.5f\t');
    end
    format=strcat(format,'\r\n');
    nBreaths=max(table.nHalfBreaths/2,1);
    data=[...
        (1:nBreaths)',...
        table.DeeFRC*100,...
        table.DeeLCI*100,...
        table.TO,...
        table.FRCAO*1e3,...
        table.DeeNormFRC*100,...
        table.DeeNormLCI*100,...
        [table.VolInsp(2:end);0]*1e3,...    % to accommodate to wBreath convention that inspiration is after expiration
        table.VolExp*1e3,...
        table.CEV*1e3,...
        table.CEV_DS*1e3,...
        table.tracerVol*1e6,...
        table.tracerVolNetto*1e6,...
        table.cumTracerVolNetto*1e6,...
        table.tracerVolReInsp*1e6,...
        table.FlowInspMean*1e6,...
        table.FlowExpMean*1e6,...
        ];
    fprintf(fid,header);
    fprintf(fid,format,data');
    
    fclose(fid);
    
    ok=1;
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting breath table related information for MBW
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag~=0
        
        maxCEV=max(table.CEV_DS);
        
        if lciStandard
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Classical method
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            infoText='  (Classical evaluation: DeeCrit = ';
            for i = 1:length(table.criticalEndRatios)
                infoText=sprintf('%s %g',infoText,table.criticalEndRatios(i)*100);
            end
            infoText=sprintf('%s %%)',infoText);
            
            cev=Direct.LCI.*Direct.FRC;
            
            hfig=figure(70+flag);
            set(hfig,'Name','LCI vs different critical Dee values');
            subplot(2,1,1);
            plot(table.CEV_DS*1e3,table.FRCAO*1e3,'+');
            maxFRCAO=max(table.FRCAO)*1e3;
            maxFRCAO=ceil(maxFRCAO*10)/10;
            hold on
            for i=1:length(cev)
                plot([cev(i)*1e3,cev(i)*1e3],[0,maxFRCAO],'k-.');
            end
            hold off
            title(strcat('FRC vs. CEV', infoText))
            xlabel('CEV / l')
            ylabel('FRC / l')
            axis([0,maxCEV*1e3,-inf,inf])
            
            subplot(2,1,2);
            semilogy(table.CEV_DS*1e3,table.DeeNormLCI*100,'+',table.CEVUpper*1e3,table.DeeNormFitUpper*100,'r-.',table.CEVLower*1e3,abs(table.DeeNormFitLower)*100,'r-.');
            hold on
            for i=1:length(cev)
                semilogy([cev(i)*1e3,cev(i)*1e3],[0.1,100],'k-.');
            end
            hold off
            title(strcat('Dee(norm) vs CEV', infoText))
            xlabel('t / s')
            xlabel('CEV / l')
            ylabel('Dee(norm) / %')
            axis([0,maxCEV*1e3,0.1,100])
        end
        
        if TOevaluation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Turnover based method
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            infoText='  (Turnover based evaluation TOcrit = ';
            for i = 1:length(table.criticalTurnovers)
                infoText=sprintf('%s %g',infoText,table.criticalTurnovers(i));
            end
            infoText=sprintf('%s)',infoText);
            
            hfig=figure(80+flag);
            set(hfig,'Name','Evaluation of Dee as function of critical TO values');
            subplot(2,1,1);
            plot(table.CEV_DS*1e3,table.FRCAO*1e3,'b+');
            maxFRCAO=max(table.FRCAO)*1e3;
            maxFRCAO=ceil(maxFRCAO*10)/10;
            title(strcat('FRC vs. TO', infoText))
            xlabel('CEV / l')
            ylabel('FRC / l')
            axis([0,maxCEV*1e3,0,maxFRCAO])
            
            subplot(2,1,2);
            semilogy(table.CEV_DS*1e3,table.DeeNormLCI*100,'+',table.CEVUpper*1e3,table.DeeNormFitUpper*100,'r-.',table.CEVLower*1e3,table.DeeNormFitLower(floor(end/2):end)*100,'r-.');
            hold on
            for i=1:length(table.Turnover.DeeNorm)
                Dee=table.Turnover.DeeNorm(i);
                semilogy([0,table.criticalTurnovers(i)*maxFRCAO],[Dee,Dee]*100,'k-.');
            end
            hold off
            title(strcat('Dee(norm) vs CEV', infoText))
            xlabel('t / s')
            xlabel('CEV / l')
            ylabel('Dee(norm) / %')
            axis([0,maxCEV*1e3,1,100])
        end
        
        if lciByFit
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fit-based method for LCI and FRC, etc. determination
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            infoText='  (critical ratios = ';
            for i = 1:length(table.criticalEndRatios)
                infoText=sprintf('%s %g',infoText,table.criticalEndRatios(i)*100);
            end
            infoText=sprintf('%s %%)',infoText);
            
            Fit=table.Fit;
            cev=Fit.LCI.*Fit.FRC;
            
            hfig=figure(90+flag);
            set(hfig,'Name','LCI based on fitted Dee curves (different critical Dee)');
            subplot(2,1,1);
            plot(table.CEV_DS*1e3,table.FRCAO*1e3,'b+');
            maxFRCAO=max(table.FRCAO)*1e3;
            maxFRCAO=ceil(maxFRCAO*10)/10;
            hold on
            for i=1:length(cev)
                if cev(i)>0
                    handleFRC=Fit.handleFRC{i};
                    cevData=0:cev(i)/100:cev(i);
                    plot(cevData*1e3,handleFRC(cevData)*1e3,'r-.');
                    plot([cev(i)*1e3,cev(i)*1e3],[0,maxFRCAO],'k-.');
                end
            end
            hold off
            title(strcat('FRC vs. CEV', infoText))
            xlabel('CEV / l')
            ylabel('FRC / l')
            axis([0,maxCEV*1e3,-inf,inf])
            
            subplot(2,1,2);
            semilogy(table.CEV_DS*1e3,table.DeeNormLCI*100,'b+');
            hold on
            for i=1:length(cev)
                if cev(i)>0
                    handleDee=Fit.handleDee{i};
                    cevData=0:cev(i)/100:cev(i)*1.5;
                    semilogy(cevData*1e3,handleDee(cevData')*100,'r-.');
                    semilogy([cev(i)*1e3,cev(i)*1e3],[1,100],'k-.');
                end
            end
            hold off
            title(strcat('Dee(norm) fit vs CEV', infoText))
            xlabel('CEV / l')
            ylabel('Dee(norm) / %')
            axis([0,maxCEV*1e3,1,100])
        end
    end
end % writeBreathTableWbreath


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Non-linear fit with 3 or 5 DOF, depending on the nature of the function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters,handle,xFit,yFit,rms,r2]=getFit(type,xMin,xMax,xInput,yInput)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Reduction of the data to the requested interval
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indices=find(1==(xInput>xMin).*(xInput<xMax));
    x=xInput(indices);
    y=yInput(indices);
    
    if strcmp(type,'nonlin')
        %if 1
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fit with novel method based on 
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [parameters,rms,r2,handle]=multiExponentialFit(x,y);
%         else %is useless if "if 1" is the other option
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % The non-linear fit is a double exponential ansatz function that is found by fitting
%             % first to a single exponential with constant term, used as initialization.
%             % Stability is satisfacorily.
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             yLog=log(y);
%             [p0,rms,r2]=linearFit(x,yLog);
%             p0=[0,exp(p0(1)),-p0(2)]';
%             singleExp=@(v,p) p(1)+p(2)*exp(-v*p(3));
%             scale=[1,1,100]'*1e-6;
%             [parametersSingle,rmsSingle,r2Single]=exponentialFit(x,y,singleExp,scale,p0);
%             handleSingle=@(v) singleExp(v,parametersSingle);
%             
%             doubleExp=@(v,p) p(1)*exp(-v*p(2))+p(3)*exp(-v*p(4));
%             scale=[1,100,1,100]'*1e-6;
%             p0=[parametersSingle(2),parametersSingle(3),parametersSingle(1),0]';
%             [parameters,rms,r2]=exponentialFit(x,y,doubleExp,scale,p0);
%             handle=@(v) doubleExp(v,parameters);
%             
%             if parameters(2)<0 || parameters(4)<0
%                 disp('WARNING: resorting to single exponential curve')
%                 parameters=parametersSingle;
%                 handle=handleSingle;
%                 rms=rmsSingle;
%                 r2=r2Single;
%             end
%         end
    else
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Conventional least square fit to linear function
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [parameters,rms,r2]=linearFit(x,y);
        handle = @(v) parameters(1)+parameters(2)*v;
    end
    
    %     fit=handle(xInput);
    xFit=x;
    yFit=handle(xFit);
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Linear fitting procedure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [parameters,rms,r2]=linearFit(x,y)

    if isempty(x) && isempty(y)
        disp('WARNING (linearFit): cannot calculate fit')
        parameters=zeros(1,2);
        rms=0;
        r2=0;
        return;
    end
    
    matrixA=[ones(size(x)),x];
    parameters=matrixA\y;
    fit=matrixA*parameters;
    
    if isempty(y)
        error('(linearFit): cannot calculate mean of y');
    else
        yMean=mean(y);
    end
    residuals=y-fit;
    variation=y-yMean;
    
    if (variation'*variation)==0
        r2=0;
    else
        r2=1-(residuals'*residuals)/(variation'*variation);
    end
    
    if isempty(y) || yMean==0
        rms=0;
        error('(linearFit): cannot calculate rms');
    else
        rms=sqrt((residuals'*residuals))/length(y)/abs(yMean);
    end
end




% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % % Find index in list
% % % % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % % % function [index]=getLCIIndex(indices)
% % % % % 
% % % % %     index=0;
% % % % %     notFound=1;
% % % % %     i=1;
% % % % %     while (i<=length(indices)-2 && notFound)
% % % % %         if ((indices(i+1)==indices(i)+1) && (indices(i+2)==indices(i)+2))
% % % % %             index=indices(i);
% % % % %             return;
% % % % %         else
% % % % %             i=i+1;
% % % % %         end
% % % % %     end
% % % % % end
% % % % % 
% % % % % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find index in list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [m0,m1,m2]=moments(TO,Dee0,TOcrit)

    targetIndices=find(TO(:,1)<=TOcrit);
    usesAll=length(targetIndices)==length(TO);
    
    if ~usesAll
        deltaTOtotal=TO(targetIndices(end)+1,1)-TO(targetIndices(end),1);
        deltaTOleft=TOcrit-TO(targetIndices(end),1);
        deltaTOright=TO(targetIndices(end)+1,1)-TOcrit;
    end
    
    m0=trapz(TO(targetIndices,1),Dee0(targetIndices,1));
    if ~usesAll
        m0=m0+(deltaTOleft*Dee0(targetIndices(end)+1,1)+deltaTOright*Dee0(targetIndices(end),1))/deltaTOtotal;
    end
    
    Dee1=Dee0.*TO;
    m1=trapz(TO(targetIndices,1),Dee1(targetIndices,1));
    if ~usesAll
        m1=m1+(deltaTOleft*Dee1(targetIndices(end)+1,1)+deltaTOright*Dee1(targetIndices(end),1))/deltaTOtotal;
    end
    
    Dee2=Dee1.*TO;
    m2=trapz(TO(targetIndices,1),Dee2(targetIndices,1));
    if ~usesAll
        m2=m2+(deltaTOleft*Dee2(targetIndices(end)+1,1)+deltaTOright*Dee2(targetIndices(end),1))/deltaTOtotal;
    end
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Evaluated tidal indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ok,table] = tidalAnalysis(table,parameters)

    nBreathMax = length(table.MMei);
    breath0 = parameters.Simulation.fromBreathTidal;
    nBreath = parameters.Simulation.nBreathTidal;
    if breath0 > 0
        startBreath = breath0;
        endBreath = breath0+nBreath-1;
    else
        endBreath = nBreathMax+breath0+1;
        startBreath = nBreathMax+breath0+2-nBreath;
    end
    if endBreath>nBreathMax || startBreath<1
        fprintf('WARNING: cannot evaluated tidal means between %d and %d, not enough breaths\n',startBreath,endBreath)
        
        table.tidalMeans.breath0            =   0;  % dummy values
        table.tidalMeans.nBreath            =   0;
        table.tidalMeans.startBreath        =   0;
        table.tidalMeans.endBreath          =   0;
        table.tidalMeans.breathRate         =   0;
        table.tidalMeans.minuteVentilation  =   0;
        table.tidalMeans.inspDuration       =   0;
        table.tidalMeans.expDuration        =   0;
        table.tidalMeans.timeInspPeakFlow   =   0;
        table.tidalMeans.timeExpPeakFlow    =   0;
        table.tidalMeans.VolTidal           =   0;
        table.tidalMeans.tPTEF2tERatio      =   0;
        ok = 0;
    else
        indicesMean = (startBreath:1:endBreath);
        table.tidalMeans.breath0            =   breath0;
        table.tidalMeans.nBreath            =   nBreath;
        table.tidalMeans.startBreath        =   startBreath;
        table.tidalMeans.endBreath          =   endBreath;
        table.tidalMeans.breathRate         =   length(indicesMean)/sum(table.breathDuration(indicesMean));
        table.tidalMeans.minuteVentilation  =   sum(table.VolExp(indicesMean))/sum(table.breathDuration(indicesMean));
        table.tidalMeans.inspDuration       =   mean(table.inspDuration(indicesMean));
        table.tidalMeans.expDuration        =   mean(table.expDuration(indicesMean));
        table.tidalMeans.timeInspPeakFlow   =   mean(table.timeInspPeakFlow(indicesMean));
        table.tidalMeans.timeExpPeakFlow    =   mean(table.timeExpPeakFlow(indicesMean));
        table.tidalMeans.VolTidal           =   mean(table.VolTidal(indicesMean));
        table.tidalMeans.tPTEF2tERatio      =   mean(table.tPTEF2tERatio(indicesMean));
        
        fprintf('\nTidal means between for breaths %d .. %d\n',startBreath,endBreath);
        fprintf('breathRateMean   =%7.2f 1/min\n',table.tidalMeans.breathRate*60);
        fprintf('minuteVentilation=%7.2f ml/min\n',table.tidalMeans.minuteVentilation*1e6*60);
        fprintf('inspDuration     =%7.2f s\n',table.tidalMeans.inspDuration);
        fprintf('expDuration      =%7.2f s\n',table.tidalMeans.expDuration);
        fprintf('timeInspPeakFlow =%7.2f s\n',table.tidalMeans.timeInspPeakFlow);
        fprintf('timeExpPeakFlow  =%7.2f s\n',table.tidalMeans.timeExpPeakFlow);
        fprintf('VolTidal         =%7.2f ml\n',table.tidalMeans.VolTidal*1e6);
        fprintf('tPTEF2tERatio    =%7.2f \n',table.tidalMeans.tPTEF2tERatio);
        ok = 1;
    end
end


