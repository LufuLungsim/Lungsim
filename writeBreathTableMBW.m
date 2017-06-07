%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Function to write ASCII breath table data
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.2, 08. Jan. 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [table] = writeBreathTableMBW(table,file,parameters,flag, signalCorrected, gasesByName, O2Washin, TracerGasWashin)

    compatibility   =   parameters.Simulation.compatibility;	% Tidal analysis
    tidal           =   parameters.Simulation.tidalAnalysis;	% Tidal analysis
    n2mbw           =   parameters.Simulation.n2mbwAnalysis;	% N2-MBW
    dtgsbw          =   parameters.Simulation.dtgsbwAnalysis;   % DTG-SBW
    dtgmbw          =   parameters.Simulation.dtgmbwAnalysis;   % DTG-MBW
    dtgmbwBaby      =   parameters.Simulation.dtgmbwBaby;
    bFile           =   parameters.Simulation.bFile;
    fastSlow        =   parameters.Simulation.fastSlow;
    torontoFile     =   parameters.Simulation.torontoFile;      % mass spectrometer based data
    MissPlexy       =   parameters.Simulation.MissPlexy;
    indexCrit       =   parameters.Operation.nVentilation;      % #breaths for minute ventilation evaluation
    expFit          =   parameters.Simulation.ExponentialByFit;
    graphState      =   parameters.Simulation.graphState;
    
    ScondSacin      =   parameters.Simulation.ScondSacin;       % Scond/Sacin evaluation
    capnoIndices    =   parameters.Simulation.CapnoIndices;
    tidalMean       =   parameters.Simulation.tidalMeanState;
    TOevaluation    =   parameters.Simulation.TOevaluation;
    lciStandard    	=   parameters.Simulation.lciStandard;
    lciByFit      	=   parameters.Simulation.lciByFit;
    
    DSpost          =   parameters.Device.volumePostcap;   
    DSpre           =   parameters.Device.volumePrecap;

    
    DS              =   DSpost+DSpre;
    nGases          =   length(gasesByName);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation of MBW breath table: determination of start cet value and
    % cropping table data to MBW breaths only (defined by the first gas
    % species)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if n2mbw || dtgmbw || fastSlow
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Critical end tidal values and turnovers
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if parameters.Simulation.fullOutput %full output selected
        table.criticalEndRatios=[2.5,3,4,5,6,7,9,12,15,18]/100;   % critical end values
        table.criticalTurnovers=[5,6,7,8];                        % critical turnover values
        end
        
        if parameters.Simulation.reducedOutput    %reduced output selected      
            table.criticalEndRatios=[2.5,5]/100;   % critical end values
            table.criticalTurnovers=[5,6,7,8]; 
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Detection of MBW begin for each gas species in the list
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1 : nGases
            [mbwBreaths, numberOfBreaths] = getMbwBreaths(table,parameters, O2Washin, TracerGasWashin, gasesByName{j});
            table.nHalfBreaths=2*length(mbwBreaths);
            [table.(gasesByName{j}).cetStart, ...
             table.(gasesByName{j}).normFactor]   =getCetStart(mbwBreaths(1),table.(gasesByName{j}),parameters, O2Washin, TracerGasWashin);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Save data for the check prephase Project returnToBaseline
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        checkPrephase = false;
        if checkPrephase
            firstMBWBreath = find(table.N2.eti<0.1, 1);
            if isempty(firstMBWBreath)
                lastPrephaseBreath = length(table.N2.et);
            else
                lastPrephaseBreath = firstMBWBreath; %No offset becasue ofo diff function
            end
            prePhaseBreaths = 1:lastPrephaseBreath;
            saveCheckPrephaseData(table, prePhaseBreaths, parameters )
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Cropping table to MBW breaths only
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        table=cropMBWBreaths(table, mbwBreaths, numberOfBreaths);
    end
            
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculation diverse diagnostic indices for each gas species
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for j = 1 : nGases
        gasByName=gasesByName{j}; 
        if n2mbw || dtgmbw || fastSlow
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set choosen gas species to variable Gas, evaluate cet_start and
            % normalized values, critical et values and some extra fields for
            % FRC, LCI and TO evluation
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gas=table.(gasByName);
            gas.General=getExtraFields(table,gas,DS,DSpre, O2Washin, TracerGasWashin);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LCI = LCI(Cet_norm) baaaased on standard criterium
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if lciStandard || fastSlow
                gas = getLCIStandard(table,gas,compatibility, parameters, TracerGasWashin);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % LCI = LCI(Cet_norm) determination using double exponential fits
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %Exclude breathvalues for fit where validData = 0
            validData=ones(size(gas.General.cev_ds));
            %Automatic outlier detection
            if parameters.Simulation.LCIbyFitAutoVal
                validData = autoDetectOutlier(gas.General.cev_ds,gas.General.cet_norm,validData,parameters);
            end
            %Manual outlier detection
            if parameters.Simulation.LCIbyFitManVal && (dtgmbw || lciByFit || fastSlow)
                validData  = manualOutlierDetection(gas.General.cev_ds,gas.General.cet_norm,validData,signalCorrected,gasByName,parameters);
            end
            
            if lciByFit %|| dtgmbw 
                gas = getFRCandLCIbyFit(validData,table,gas,compatibility, parameters, TracerGasWashin);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Turnover based detection of Cet(norm)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if TOevaluation
                gas = getEtFromTurnover(table,gas);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fitting linear slopes in semilog plot (used only as a guide for the eye in graphs
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gas.General.logCet_norm = log(gas.General.cet_norm);
            CEVmax = max(gas.General.cev_ds);
            [~,~,fitUpper,~,~]=getFit('lin',0,         CEVmax/3,gas.General.cev_ds,gas.General.logCet_norm);
            [~,~,fitLower,~,~]=getFit('lin',2*CEVmax/3,CEVmax,  gas.General.cev_ds,gas.General.logCet_norm);
            gas.General.cet_normFitUpper = exp(fitUpper);
            gas.General.cet_normFitLower = exp(fitLower);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Evaluation of Scond/Sacin and data to plot the corresponding figure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if ScondSacin
                lBound  =   parameters.Simulation.lBound;
                uBound  =   parameters.Simulation.uBound;
                gas     =   evaluateScondSacin(1,lBound,uBound,gas,table,parameters);
                uBound	=   parameters.Simulation.uBoundStar;
                lBound	=   parameters.Simulation.lBoundStar;
                gas     =   evaluateScondSacin(2,lBound,uBound,gas,table,parameters);
            end
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fast Slow Analysis
            % First decompose the signal into fast and slow ratios by
            % fitting the signal to mathematical function.
            % Secondly find the assumed volume of the fast and slow
            % compartments as well as the diffusion factor between the two
            % by fitting the measured signal to a simplyfied model of the
            % lung
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if fastSlow 
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                gas = calculateFastSlow(table, gas, validData, parameters);
                %gas = fitToLungModel(table, gas, parameters);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Turnover based detection of Cet(norm) using the double exponential fit
            % Attention: this part is omitted, as the TO is known only implicitly within this setup
            % It is possible to solve this issue, but due to limited resources, this will be done later.
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
        elseif tidal
            gas = table.(gasByName);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Set dummy values (to facilitate postprocessing)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            gas.General.cet_crit=zeros(size(gas.et));
            gas.General.to=zeros(size(gas.et));
            gas.General.frcao=zeros(size(gas.et));
            gas.General.cet_start=zeros(size(gas.et));
            gas.General.cet_norm=zeros(size(gas.et));
            gas.General.cev=zeros(size(gas.et));
            gas.General.cev_ds=zeros(size(gas.et));
            gas.insp=zeros(size(gas.et));
            gas.exp=zeros(size(gas.et));
            gas.General.netExpVol=zeros(size(gas.et));
            gas.General.cumNetExpVol=zeros(size(gas.et));
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Calculating fields for tidal breathing maneuvers
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            if capnoIndices
                table.Capno = getCapnoIndices(table,parameters);
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Write breath table
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        writeBreathTableData(file, table, gas, parameters);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Reset Gas record
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        table.(gasByName)=gas;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of capno indices (only once, not for each gas)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if capnoIndices
        table.Capno = getCapnoIndices(table,parameters);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of tidal averages, etc.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    table.Tidal = gedTidalInfo(table, indexCrit);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Calculate tidal means
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (tidal || n2mbw || dtgmbw ) && tidalMean
        table.TidalMeans = getTidalMeans(table, parameters);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting breath table related information for dtgsbw and n2mbw
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag~=0 && dtgsbw
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot charaterizing the sbw
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        hfig=figure(40+flag);
        set(hfig,'Name','SBW characterization');
        
        figureSize=get(hfig,'Position');    % rescaling graph if necessary (for better visibility)
        height=figureSize(4);
        width=figureSize(3);
        if width>height
            figureSize(4)=1.6*height;
            figureSize(2)=figureSize(2)-0.6*height;
            set(hfig,'Position',figureSize);
        end
        
        maxVol=max([table.vExp;table.vInsp])*1e6;
        maxVolExp=max(table.vExp)*1e6;
        volIntercept=table.InterceptionVolume(1)*maxVolExp;
        volDiffMax=table.DiffMMssMaxVolume(1)*maxVolExp;
        
        subplot(3,1,1);
        plot(table.vExp*1e6,table.MMss*1e3,'b');
        title('MMss vs Vol_e_x_p')
        xlabel('Vol / ml')
        ylabel('MMss / g/mol')
        axis([0,maxVol,-inf,inf])
        
        maxDiffMMss=max(abs(table.DiffMMss))*1e3;
        diffMMssMax=table.DiffMMssMax(1)*1e3;
        subplot(3,1,2);
        plot(table.vExp*1e6,table.DiffMMss*1e3,'b');
        hold on
        plot([volIntercept,volIntercept],[-maxDiffMMss/5,maxDiffMMss/5],'-.k');
        plot([volIntercept-maxVol/20,volIntercept+maxVol/20],[0,0],'-.k');
        text(volIntercept*1.05,0.04,sprintf('vInter=%2.1f%%',table.InterceptionVolume(1)*100) )
        plot([volDiffMax,volDiffMax],[diffMMssMax-maxDiffMMss/5,diffMMssMax+maxDiffMMss/5],'-.k');
        plot([volDiffMax-maxVol/20,volDiffMax+maxVol/20],[diffMMssMax,diffMMssMax],'-.k');
        text(volDiffMax*1.05,diffMMssMax-0.08,sprintf('vMax=%2.1f%%\ndiffMax=%.2f',table.DiffMMssMaxVolume(1)*100,diffMMssMax) )
        hold off
        title('DiffMMss vs Vol_e_x_p')
        xlabel('Vol / ml')
        ylabel('DeltaMMss / g/mol')
        axis([0,maxVol,-0.4,0.4])
        
        subplot(3,1,3);
        plot(table.vExp*1e6,table.IvExp*1e6,'b');
        hold on
        plot((table.vInsp(end)-table.vInsp)*1e6,table.IvInsp*1e6,'r');
        hold off
        title(strcat('Volume flow diagram for SBW'));
        xlabel('vol / ml')
        ylabel('Iv / ms/s')
        legend('exp','insp')
        axis([0,maxVol,-inf,inf])
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MBW output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if flag~=0 && n2mbw || dtgmbw
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot capno indices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if parameters.Simulation.CapnoIndices
            plotCapnoIndices(table,table.Capno,[0 0 1],100+flag);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot diagnostic information by gas species
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for j = 1 : nGases
            gasByName=gasesByName{j};
            gas = table.(gasByName);
            if graphState
                plotScondSacin(1,table,gas,parameters,[0 1 0],j*1000+110+flag);
                plotScondSacin(2,table,gas,parameters,[1 0 1],j*1000+120+flag);
            end
            
            maxCEV=max(gas.General.cev_ds);
            
            if lciStandard
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Classic method
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                infoText=sprintf('Classic LCI for %s, C_e_t_,_C_r_i_t=',gasByName);
                for i = 1:length(table.criticalEndRatios)
                    infoText=sprintf('%s %g',infoText,table.criticalEndRatios(i)*100);
                end
                infoText=sprintf('%s %%',infoText);
                
                cev=gas.Standard.lci.*gas.Standard.frc;
                
                if graphState
                    hfig=getOrMakeFigure('LCI vs different critical Cet values');
                    subplot(2,1,1);
                    plot(gas.General.cev_ds*1e3,gas.General.frcao*1e3,'+');
                    maxFRCAO=max(gas.General.frcao)*1e3;
                    maxFRCAO=ceil(maxFRCAO*10)/10;
                    hold on
                    for i=1:length(cev)
                        plot([cev(i)*1e3,cev(i)*1e3],[0,maxFRCAO],'k-.');
                    end
                    hold off
                    title(infoText)
                    xlabel('CEV / l')
                    ylabel('FRC / l')
                    axis([0,maxCEV*1e3, -inf,inf])

                    warning('off','MATLAB:Axes:NegativeDataInLogAxis')
                    warning('off','MATLAB:plot:IgnoreImaginaryXYPart')
                    subplot(2,1,2);
                    semilogy(gas.General.cev_ds*1e3,gas.General.cet_norm*100,'+',gas.General.cev_ds*1e3,gas.General.cet_normFitUpper*100,'r-.',gas.General.cev_ds(floor(end/2):end)*1e3,gas.General.cet_normFitLower(floor(end/2):end)*100,'r-.');
                    hold on
                    for i=1:length(cev)
                        semilogy([cev(i)*1e3,cev(i)*1e3],[1,100],'k-.');
                    end
                    %2.5%marker line
                    line('XData', [0 maxCEV*1e3], 'YData', [2.5 2.5], 'LineStyle', '--', 'Color', 'g');
                    hold off
                    xlabel('CEV / l')
                    ylabel('Cet(norm) / %')
                    axis([0,maxCEV*1e3,1,100])
                    warning('on','MATLAB:Axes:NegativeDataInLogAxis')
                    warning('on','MATLAB:plot:IgnoreImaginaryXYPart')
                end
            end
            
            if TOevaluation
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Turnover based method
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                infoText='Turnover based evaluation TOcrit = ';
                for i = 1:length(table.criticalTurnovers)
                    infoText=sprintf('%s %g',infoText,table.criticalTurnovers(i));
                end
                infoText=sprintf('%s',infoText);
                
                hfig=figure(j*1000+80+flag);
                set(hfig,'Name','Evaluation of Cet as function of critical TO values');
                subplot(2,1,1);
                plot(gas.General.cev_ds*1e3,gas.General.frcao*1e3,'b+');
                maxFRCAO=max(gas.General.frcao)*1e3;
                maxFRCAO=ceil(maxFRCAO*10)/10;
                title(infoText)
                xlabel('CEV / l')
                ylabel('FRC / l')
                axis([0,maxCEV*1e3,0,maxFRCAO])
                
                subplot(2,1,2);
                semilogy(gas.General.cev_ds*1e3,gas.General.cet_norm*100,'+',gas.General.cev_ds*1e3,gas.General.cet_normFitUpper*100,'r-.',gas.General.cev_ds(floor(end/2):end)*1e3,gas.General.cet_normFitLower(floor(end/2):end)*100,'r-.');
                hold on
                for i=1:length(gas.Turnover.cet_norm)
                    cet=gas.Turnover.cet_norm(i);
                    semilogy([0,table.criticalTurnovers(i)*maxFRCAO],[cet,cet]*100,'k-.');
                end
                %2.5%marker line
                line('XData', [0 maxCEV*1e3], 'YData', [2.5 2.5], 'LineStyle', '--', 'Color', 'g');
                hold off
                title(strcat('Cet(norm) vs CEV', infoText))
                xlabel('t / s')
                xlabel('CEV / l')
                ylabel('Cet(norm) / %')
                axis([0,maxCEV*1e3,1,100])
            end
            
            if lciByFit
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Fit-based method for LCI and FRC, etc. determination
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                infoText=sprintf('By fit LCI for %s, C_e_t_,_C_r_i_t= ',gasByName);
                for i = 1:length(table.criticalEndRatios)
                    infoText=sprintf('%s %g',infoText,table.criticalEndRatios(i)*100);
                end
                infoText=sprintf('%s %%',infoText);
                
                cev=gas.ByFit.lci.*gas.ByFit.frc;
                
                hfig=figure(j*1000+90+flag);
                set(hfig,'Name','LCI by fit of Cet curves (different critical Cet)');
                
                figureSize=get(hfig,'Position');    % rescaling graph if necessary (for better visibility)
                height=figureSize(4);
                width=figureSize(3);
                if width>height
                    figureSize(4)=2.1*height;
                    figureSize(2)=figureSize(2)-1.1*height;
                    set(hfig,'Position',figureSize);
                end

                warning('off','MATLAB:Axes:NegativeDataInLogAxis')
                warning('off','MATLAB:plot:IgnoreImaginaryXYPart')
                subplot(3,1,1);
                semilogy(gas.General.cev_ds*1e3,gas.General.cet_norm*100,'b+');
                hold on
                for i=1:length(cev)
                    if cev(i)>0
                        k=gas.ByFit.nBreaths(i);
                        n=find(validData(1:k-1)==1, 1, 'last');
                        if expFit
                            semilogy(gas.General.cev_ds*1e3,gas.ByFit.handleCet{1}(gas.General.cev_ds)*1e2,'k--')
                        else
                            semilogy([gas.General.cev_ds(n)*1e3, gas.General.cev_ds(k)*1e3], [gas.General.cet_norm(n)*100, gas.General.cet_norm(k)*100],'k--');
                        end
                        semilogy([cev(i)*1e3,cev(i)*1e3],[1,100],'k-.');
                    end
                end
                %2.5%marker line
                line('XData', [0 maxCEV*1e3], 'YData', [2.5 2.5], 'LineStyle', '--', 'Color', 'g');
                hold off
                title(infoText)
                xlabel('CEV / l')
                ylabel('Cet(norm) / %')
                axis([0,maxCEV*1e3,1,100])
                warning('on','MATLAB:Axes:NegativeDataInLogAxis')
                warning('on','MATLAB:plot:IgnoreImaginaryXYPart')
                
                hSubplot30(1)=subplot(3,1,2);
                if torontoFile
                    plot(signalCorrected.ts,signalCorrected.CO2*100, signalCorrected.ts,signalCorrected.He*100, signalCorrected.ts,signalCorrected.SF6*100,'Parent',hSubplot30(1));
                    hold on
                    %plot Breaths
                    breathTimes=signalCorrected.ts(signalCorrected.breathIndex);
                    for i=1:length(breathTimes)
                        plot([breathTimes(i), breathTimes(i)],[0,100],'k-.');
                    end
                    hold off
                    title('molar fractions at POI')
                    xlabel('t / s')
                    ylabel('x / %')
                    legend('CO2', 'He', 'SF6')
                    axis([-inf,inf,-inf,10])
                else
                    plot(signalCorrected.ts,signalCorrected.CO2Poi1*100, signalCorrected.ts,signalCorrected.O2Poi1*100, signalCorrected.ts,signalCorrected.N2Poi1*100,...
                        signalCorrected.ts,signalCorrected.ArPoi1*100,'Parent',hSubplot30(1));
                    hold on
                    %plot Breaths
                    breathTimes=signalCorrected.ts(signalCorrected.breathIndex);
                    
                    for i=1:length(breathTimes)
                        plot([breathTimes(i), breathTimes(i)],[0,100],'k-.');
                    end
                   
                    hold off
                    title('molar fractions at POI1')
                    xlabel('t / s')
                    ylabel('x / %')
                    legend('CO2', 'O2', 'N2', 'Ar')
                    axis([-inf,inf,-inf,inf])
                end
                
                hSubplot30(2)=subplot(3,1,3);
                plot(signalCorrected.ts,signalCorrected.IvEffBTPS*1e6,signalCorrected.ts,signalCorrected.VolBTEffBTPS*1e6,'Parent',hSubplot30(2));
                hold on
                for i=1:length(breathTimes)
                    plot([breathTimes(i), breathTimes(i)],[-500,500],'k-.');
                end
                hold off
                title('flow related data')
                xlabel('t / s')
                ylabel('I_V_,_B_T_P_S, V_B_T_P_S / ml/s, ml')
                legend('IvEffBTPS', 'VolBTEffBTPS')
                axis([-inf,inf,-inf,inf])
                linkaxes(hSubplot30,'x')
                
            end
            
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Tidal breaths
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    elseif flag~=0 && tidal
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plot capno indices
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if parameters.Simulation.CapnoIndices
            plotCapnoIndices(table,table.Capno,[0 0 1],100+flag);
        end
        
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Detection of MBW begin
% 
% The Washout parameter is used to detect if the signal contains a Washin
% and Washout measurement in which case the breaths have to be segmented in
% a different way. 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mbwBreaths, numberOfBreaths] = getMbwBreaths(table,parameters, O2Washin, TracerGasWashin, gasName)

    dtgmbwBaby         =   parameters.Simulation.dtgmbwBaby;     %SF6 4% %%TODO rename dtgmbwBaby to dtgmbw4SF6

    if ~parameters.Simulation.torontoFile                     % spiroware based data set
        %Find the beginning of O2 washin which should coincide with a
        %Molarmass increase
        [~, StartMMWashin]    = max(diff(table.MMss.eti));     %Find the biggest increase in MMss
        [~, StartMMWashout]   = min(diff(table.MMss.eti));     %Find the biggest decrease in MMss
           
        %If the prephase is not stable probably the file was cut
        %incorrectly
        etFirst3Breaths = table.MMss.et(1:3)./mean(table.MMss.et(1:3));
        X = [ones(length(etFirst3Breaths),1) [1:length(etFirst3Breaths)]']\etFirst3Breaths;
        if abs(X(2)) > 0.025
            StartMMWashout = 0;
            StartMMWashin  = 0;
        end
        
        %(+1 to offset the shift caused by the diff() function)     
        StartMMWashin         = StartMMWashin +1;              %Compensate for index shift through diff()
        
        %In a normal routine O2 is Washed in and then out. Depending if the
        %tracer gas is naturally occuring (N2) or is washed in with the O2,
        %the relative washin/washout phase is set here
        if parameters.Simulation.dtgmbwAnalysis
            if	    O2Washin
                if parameters.Simulation.aFile
                    mbwBreaths = StartMMWashin:StartMMWashout;
                else 
                    mbwBreaths = StartMMWashin:length(table.MMss.et);
                end
            elseif ~O2Washin
                StartMMWashout	= StartMMWashout +1;              %Compensate for index shift through diff()
                mbwBreaths      = StartMMWashout:length(table.MMss.et);
            end
        else
            mbwBreaths = StartMMWashin:length(table.MMss.et);
        end
            
        if isempty(mbwBreaths)
            warning(['(writeBreathTableMBW): cannot determine MBW breaths, using all of the measurement for the analysis. Check the input file %s !!!', parameters.Simulation.fileInput]);
            mbwBreaths=1:length(table.MMss.insp(:,1));
        end
    else                                    % mass spectrometer based data set
        normConcentration=parameters.Operation.toronto(3);  % assuming equal start of He and SF6
        mbwBreaths=find(table.He.meanInsp(:,1)<normConcentration/2);
    end
    
    numberOfBreaths = length(table.MMss.et);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Set Cet_Start (normalization of signal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [cetStart, normFactor] = getCetStart(startIndex, gas, parameters, O2Washin, TracerGasWashin)
    
    %If there is at least 1 breath before the start of the Washin
    if startIndex>=2
        %If there is one breath before the start of the Washin set out a
        %warning
        if startIndex == 2
            warning('(writeBreathTableMBW): cannot determine Cet_start precisely, only one pre-breath found');    
        end
        %Use 1-3 breaths to determine the cetStart
        cetStart = mean(gas.et(max(startIndex-3, 1):startIndex-1,1)); 
        
        if TracerGasWashin
            normFactor = max(mean(gas.eti(startIndex+1:startIndex+6)), cetStart);
        else
            normFactor = cetStart;
        end
    %If no prephase breath was found use standard values
    else
        warning('(writeBreathTableMBW): cannot determine Cet_start precisely, no pre-breath found');
        %If the tracer gas is beeing washed in, use 0% as starting
        %concentration
        if TracerGasWashin 
            cetStart = 0;
            normFactor = max(mean(gas.eti(startIndex+1:startIndex+6)), max(gas.et));
        %Else use natural constants or normed gas values
        else
            if strcmp(gas.name, 'N2')
                cetStart = parameters.Operation.xAir(5);
            elseif strcmp(gas.name, 'SF6')
                if parameters.Simulation.dtgmbwBaby;
                    cetStart = 0.04;
                else
                    cetStart = 0.03;
                end    
            else
                error(['(writeBreathTableMBW): cannot determine Cet_start precisely, Gas not recognized' gas.name])
            end
            normFactor = cetStart;
        end
    end
    cetStart = max(cetStart, 0);   
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating extra fields, FRC, LCI, etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [general] = getExtraFields(table, gas, DS, DSpre, O2Washin, TracerGasWashin)
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % MBW related information related evaluations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    cetStart            = gas.cetStart;
    general.cet_start   = ones(size(gas.et))*cetStart;
    
    if TracerGasWashin
        general.cet_norm     = 1-gas.et./gas.normFactor;
        general.netExpVol    = gas.expRoof-gas.inspRoof;     
    else
        general.cet_norm    = gas.et./gas.normFactor;
        general.netExpVol   = gas.exp-gas.insp;                   % reInsp calculation with N2.et*DSpost        
    end
    
    general.cet_crit    = table.criticalEndRatios;	% set of actual cet_crit
    general.cev=cumsum(table.volExp); %partialvolumen
    general.cev_ds=general.cev-cumsum(ones(size(table.volExp)))*DS; 
    general.cumNetExpVol=cumsum(general.netExpVol);
    general.frcsp=general.cumNetExpVol./abs(general.cet_start-gas.et);
    general.frcao=general.frcsp-DSpre;    
    general.to=general.cev_ds./general.frcao;
    general.minuteVentilation=general.cev_ds./cumsum(table.breathDuration);
    % TODO es ist unklar, was hier als minuteVentilation definiert sein soll:
    % ist es das tatsächlich ausgeatmete Luftvolumen (cev), oder dasjenige,
    % welches um den Totraum vermindert wurde (wie es für die TO-Berechnung
    % verwendet wird)? Hier ist es so gemacht, dass es mit der
    % getFRCandLCIbyFit-Methode übereinstimmt!!
    % general.minuteVentilation=general.cev./cumsum(table.breathDuration);
    
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Cet(norm) based determination of TO=LCI
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gas] = getLCIStandard(table, gas, compatibility, parameters, Washin)

    momentRatio = parameters.Simulation.MomentRatio;

    gas.Standard.cet_norm=zeros(size(table.criticalEndRatios));
    gas.Standard.lci=zeros(size(table.criticalEndRatios));
    gas.Standard.frc=zeros(size(table.criticalEndRatios));
    gas.Standard.duration=zeros(size(table.criticalEndRatios));
    gas.Standard.minuteVentilation=zeros(size(table.criticalEndRatios));
    gas.Standard.breathRate=zeros(size(table.criticalEndRatios));
    gas.Standard.nBreaths=zeros(size(table.criticalEndRatios));

    targetIndexConcentration=zeros(size(table.criticalEndRatios));
    for i = 1:length(table.criticalEndRatios)
        targetIndices=find(gas.General.cet_norm<table.criticalEndRatios(i));
        targetIndexConcentration(i)=getLCIIndex(targetIndices,compatibility);
    end
           
    disp(' ');
    if Washin 
        mode = 'Washin';
    else
        mode = 'Washout';
    end
    disp(['FRC and LCI for given critical Cet(norm) of ' gas.name ' for ' mode]);
    for i = 1:length(table.criticalEndRatios)
        CetCrit=table.criticalEndRatios(i);
        if targetIndexConcentration(i)>0
            index=targetIndexConcentration(i);
            gas.Standard.cet_norm(i)=CetCrit;
            gas.Standard.lci(i)=gas.General.to(index,1);
            gas.Standard.frc(i)=gas.General.frcao(index,1);
            gas.Standard.duration(i)=sum(table.breathDuration(1:index,1));
            gas.Standard.minuteVentilation(i)=gas.General.minuteVentilation(index);
            gas.Standard.nBreaths(i)=index;
            gas.Standard.breathRate(i)=index/gas.Standard.duration(i);
            
            data=gas.Standard;
            fprintf('#breath=%3g: duration=%5.1f s, FRC=%6.1f ml, LCI=%5.2f at Cet(norm)=%4.1f %%, BR=%4.1f 1/min, MV=%6.1f ml/min\n', ...
                    data.nBreaths(i),data.duration(i),data.frc(i)*1e6,data.lci(i),CetCrit*100,data.breathRate(i)*60,data.minuteVentilation(i)*1e6*60);
            if momentRatio
                [m0, m1, m2] = getMomentRatios(gas.General.cet_norm(1:index), gas.General.to(1:index,1));
                gas.Standard.momentRatio(i,1:3) = [m0, m1, m2];
            end        
        else
            fprintf(['WARNING: cannot determine FRC and LCI for Cet(norm) = %f %% (' gas.name '  stays above thershold)\n'],CetCrit*100);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Turnover based detection of Cet(norm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gas] = getEtFromTurnover(table, gas)

    gas.Turnover.cet_norm=zeros(size(table.criticalTurnovers));
    gas.Turnover.lci=zeros(size(table.criticalTurnovers));
    gas.Turnover.frc=zeros(size(table.criticalTurnovers));
    gas.Turnover.duration=zeros(size(table.criticalTurnovers));
    gas.Turnover.nBreaths=zeros(size(table.criticalTurnovers));

    disp(' ');
    disp('FRC and Cet(norm) for given TO');
    for i = 1:length(table.criticalTurnovers)
        toCrit = table.criticalTurnovers(i);
        targetIndices=find(gas.General.to(:,1)>toCrit);
        if ~isempty(targetIndices)
            index = targetIndices(1);    % due to cropping, the startindex is now generally equal to 1
            gas.Turnover.cet_norm(i)=gas.General.cet_norm(index);
            gas.Turnover.lci(i)=toCrit;
            gas.Turnover.frc(i)=gas.General.frcao(index);
            gas.Turnover.duration(i)=sum(table.breathDuration(1:index,1));
            gas.Turnover.nBreaths(i)=index;
            
            data = gas.Turnover;
            fprintf('#breath=%d: duration=%.1f s, FRC=%.1f ml, Cet(norm)=%.2f %% at TO=%g\n', data.nBreaths(i),data.duration(i),data.frc(i)*1e6,data.cet_norm(i)*100,toCrit);
        else
            fprintf('WARNING: cannot determine Cet(norm) and FRC for TO=%g (TOmax<TOcrit)\n',toCrit);
        end
    end
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculating fields for tidal breathing maneuvers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [tidal] = gedTidalInfo(table, indexCrit)
    tidal.minuteVentilation=table.volExp./table.breathDuration;
    cev=cumsum(table.volExp);

    timeLine=cumsum(table.breathDuration);
    tidal.minuteVentilationTotal=cev(end,1)/timeLine(end,1);
    tidal.breathRateTotal=length(cev(:,1))/timeLine(end,1);
    if indexCrit<=length(cev(:,1));
        tidal.minuteVentilationCrit=cev(indexCrit,1)/timeLine(indexCrit,1);
        tidal.breathRateCrit=indexCrit/timeLine(indexCrit,1);
    else
        fprintf('WARNING: cannot determine minute ventilation, not enought breaths\n')
        tidal.minuteVentilationCrit=0;
        tidal.breathRateCrit=0;
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Tidal information related evaluations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tidal.volInspMean=mean(table.volInsp);
    tidal.volExpMean =mean(table.volExp );
    tidal.volInspStd =std(table.volInsp,1)/tidal.volInspMean;
    tidal.volExpStd  =std(table.volExp ,1)/tidal.volExpMean;

    fprintf('\nmean(VolInsp)=%.1f ml+-%.1f%%, mean(VolExp)=%.1f ml+-%.1f%%\n', tidal.volInspMean*1e6,tidal.volInspStd*100,tidal.volExpMean*1e6,tidal.volExpStd*100);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate tidal means and some output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function tidalMeans = getTidalMeans(table, parameters)

    nBreathMax = length(table.volExp);
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

        tidalMeans.breath0            =   0;  % dummy values
        tidalMeans.nBreath            =   0;
        tidalMeans.startBreath        =   0;
        tidalMeans.endBreath          =   0;
        tidalMeans.breathRate         =   0;
        tidalMeans.minuteVentilation  =   0;
        tidalMeans.inspDuration       =   0;
        tidalMeans.expDuration        =   0;
        tidalMeans.timeInspPeakFlow   =   0;
        tidalMeans.timeExpPeakFlow    =   0;
        tidalMeans.volTidal           =   0;
        tidalMeans.tPTEF2tERatio      =   0;
    else
        indicesMean = (startBreath:1:endBreath);
        tidalMeans.breath0            =   breath0;
        tidalMeans.nBreath            =   nBreath;
        tidalMeans.startBreath        =   startBreath;
        tidalMeans.endBreath          =   endBreath;
        tidalMeans.breathRate         =   length(indicesMean)/sum(table.breathDuration(indicesMean));
        tidalMeans.minuteVentilation  =   sum(table.volExp(indicesMean))/sum(table.breathDuration(indicesMean));
        tidalMeans.inspDuration       =   mean(table.inspDuration(indicesMean));
        tidalMeans.expDuration        =   mean(table.expDuration(indicesMean));
        tidalMeans.timeInspPeakFlow   =   mean(table.timeInspPeakFlow(indicesMean));
        tidalMeans.timeExpPeakFlow    =   mean(table.timeExpPeakFlow(indicesMean));
        tidalMeans.volTidal           =   mean(table.volTidal(indicesMean));
        tidalMeans.tPTEF2tERatio      =   mean(table.tPTEF2tERatio(indicesMean));

        fprintf('\nTidal means between for breaths %d .. %d\n',startBreath,endBreath);
        fprintf('breathRateMean   =%7.2f 1/min\n',tidalMeans.breathRate*60);
        fprintf('minuteVentilation=%7.2f ml/min\n',tidalMeans.minuteVentilation*1e6*60);
        fprintf('inspDuration     =%7.2f s\n',tidalMeans.inspDuration);
        fprintf('expDuration      =%7.2f s\n',tidalMeans.expDuration);
        fprintf('timeInspPeakFlow =%7.2f s\n',tidalMeans.timeInspPeakFlow);
        fprintf('timeExpPeakFlow  =%7.2f s\n',tidalMeans.timeExpPeakFlow);
        fprintf('VolTidal         =%7.2f ml\n',tidalMeans.volTidal*1e6);
        fprintf('tPTEF2tERatio    =%7.2f \n',tidalMeans.tPTEF2tERatio);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate capno indices and writing some output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function capno = getCapnoIndices(table, parameters)

        capno = evaluateCapnoIndices(table,parameters);
        disp(' ');
        disp('Capno indices: mean(std)');
        fprintf('FowlerDS=%7.2f(%7.2f) ml\n',  [capno.meanFowlerDead,capno.stdFowlerDead]*1e6);
        fprintf('SnIICO2 =%7.2f(%7.2f) 1/L\n', [capno.meanCoeffII,   capno.stdCoeffII]/1e3   );
        fprintf('SnIIICO2=%7.2f(%7.2f) 1/L\n', [capno.meanCoeffIII,  capno.stdCoeffIII]/1e3  );
        fprintf('FCO2E   =%7.2f(%7.2f) %%\n',  [capno.meanfCO2E,     capno.stdfCO2E]*1e2     );
        fprintf('CO2et   =%7.2f(%7.2f) %%\n',  [capno.meanCO2et,     capno.stdCO2et]*1e2     );
        
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write breath table into a file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function writeBreathTableData(file, table, gas, parameters)
    fid = fopen(file, 'w');
    if fid<1
        error('cannot open file %s\n',file);
    end
    
    dtgsbw     =parameters.Simulation.dtgsbwAnalysis;
    torontoFile=parameters.Simulation.torontoFile;
    
    header1='Breath #\tCet [%%]\tTO\tFRC [l]\tSnIIIms\tCet_start\tCet(norm)\tC mean slope\tC mean breath\tVolInsp [l]\tVolExp [l]\t'; % 11
    header2='CEV [l]\tCEV-DS [l]\tVolN2Insp [ml]\tVolN2Exp [ml]\tVolN2Netto [ml]\tCumVolN2Netto [ml]\tVol N2Reinsp [ml]\t'; % 7
    header3='SnIII\tSnIII, C breath\tVdCO2\tFlowInsp.mean [ml/s]\tFlowExp.mean [ml/s]\tRR\t'; % 6
    if dtgsbw   % slopes are not normed (i.e., in natural units)
        header4='SlopeCO2 [%%/l]\tSlopeN2 [%%/l]\tSlopeMMss [g/mol/l]\tSlopeMMssCalc [g/mol/l]\tSlopeDiffMMss [g/mol/l]\t'; % 5
        xSlopes=100/1000;    % %/l
        mmSlopes=1e3/1000;   % g/mol/l
    else        % slopes are normed
        header4='SlopeCO2 [1/l]\tSlopeN2 [1/l]\tSlopeMMss [1/l]\tSlopeMMssCalc [1/l]\tSlopeDiffMMss [1/l]\t'; % 5
        xSlopes=1/1000;      % 1/l
        mmSlopes=1/1000;     % 1/l
    end
    header5='InterceptionVolume [%%ExpVol]\tEndTidalDiffMMss [g/mol]\tDiffMMssMax [g/mol]\tDiffMMssMaxVolume [%%ExpVol]\t'; % 4
    header = strcat(header1,header2,header3,header4,header5,'\r\n');
    format='';
    nColumns=33; % = 11+7+6+5+4
    for i=1:nColumns
        format=strcat(format,'%.5f\t');
    end
    format=strcat(format,'\r\n');
    if dtgsbw
        nBreaths=max(floor((table.nHalfBreaths+1)/2),1);
    else
        nBreaths=max(table.nHalfBreaths/2,1);
    end
    
    zeroData=zeros(nBreaths,1);
    if ~torontoFile
        slopesCO2=table.CO2.slope(:,2)*xSlopes;
        slopesMMss=table.MMss.slope(:,2)*mmSlopes;
        slopesMMssCalc=table.MMss.calcSlope(:,2)*mmSlopes;
        slopesMMssDiff=table.MMss.diffSlope(:,2)*mmSlopes;
        interceptionVolume=table.InterceptionVolume(:,1)*100;
        endTidalDiffMMss=table.EndTidalDiffMMss(:,1)*1e3;
        diffMMssMax=table.DiffMMssMax(:,1)*1e3;
        diffMMssMaxVolume=table.DiffMMssMaxVolume(:,1)*100;
    else
        slopesMMss=zeroData;
        slopesMMssCalc=zeroData;
        slopesMMssDiff=zeroData;
        interceptionVolume=zeroData;
        endTidalDiffMMss=zeroData;
        diffMMssMax=zeroData;
        diffMMssMaxVolume=zeroData;
        slopesCO2=zeroData;
    end
    
    data=[...
        (1:nBreaths)',...
        gas.et*100,...
        gas.General.to,...
        gas.General.frcao*1e3,...
        zeroData,...
        gas.General.cet_start*100,...
        gas.General.cet_norm*100,...
        zeroData,...
        zeroData,...
        table.volInsp*1e3,...
        table.volExp*1e3,...
        gas.General.cev*1e3,...
        gas.General.cev_ds*1e3,...
        gas.insp*1e6,...
        gas.exp*1e6,...
        gas.General.netExpVol*1e6,...
        gas.General.cumNetExpVol*1e6,...
        gas.reInsp*1e6,...
        zeroData,...
        zeroData,...
        zeroData,...
        zeroData,...
        zeroData,...
        zeroData,...
        slopesCO2,...
        gas.slope(:,2)*xSlopes,...
        slopesMMss,...
        slopesMMssCalc,...
        slopesMMssDiff,...
        interceptionVolume,...
        endTidalDiffMMss,...
        diffMMssMax,...
        diffMMssMaxVolume,...
        ];
    fprintf(fid,header);
    fprintf(fid,format,data');
    
    fclose(fid);
end