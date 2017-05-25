%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading breathing manoeuver data (A-files and older or wbreath export .txt-files
%
% Version 3.0, 25. Juli 2012
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header,data,signal]=readSignalSpirowareMassSpec(file,figureNo,parameters)

    
   
    sampleFlow      =   parameters.Operation.sampleFlow;
    mmAir           =   parameters.Operation.MMair;
    mm2Percent      =   parameters.Operation.MM2Percent;
    
    torontoFile     =   parameters.Simulation.torontoFile;
    
    spirowareType   =   parameters.Operation.spirowareType;
    aFile           =   parameters.Simulation.aFile;
    bFile           =   parameters.Simulation.bFile;
    sbwMm           =	parameters.Simulation.sbwMm;
    sbwPercent      =	parameters.Simulation.sbwPercent;
    dtgmbw          =   parameters.Simulation.dtgmbwAnalysis;
    heMbw           =   parameters.Simulation.heMbw;
    MissPlexy       =   parameters.Simulation.MissPlexy;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read file and reformat raw data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
    [header,data]   =   readFile(file);               	% reading data file
    token           =   splitToken(header);                  % header tags
        
    iTime={'Time', 'Tid'};
    iFlow={'Flow', 'Flow_BTPS'};
    iHe={'He'};
    iSF6={'SF6'};
    iCO2={'CO2', 'CO2_', 'CO2_RAW', 'CO2*'};
    iMMss={'MMss'};
    iMMms={'MMms'};
    iMM={'MM'};
    iO2={'O2'};
    iSampleFlow={'SampleFlow'};
    iN2={'N2'};
    iAr={'Ar', 'AR'};
    iO2Delay={'DelayO2'};
    iDelayMMss={'DelayMMss'};
    tokens      = {iTime, iFlow, iHe, iSF6, iCO2, iMMss, iMMms, iMM, iO2, iSampleFlow, iN2, iAr, iO2Delay, iDelayMMss};
     
    %data = reshapeData(aFile, token, bFile, data, file, parameters);
    
    %Fill all available Gases, the ones not present in the file are filled
    %with 0
    for i=1:length(tokens)  
        index = findToken(tokens{i}, token, 1); 
        if index > 0 && index <= size(data, 2)      % If a token was found and there is data for the token (which is not always the case in DTGMBW Files)
            signal.(matlab.lang.makeValidName(tokens{i}{1}))=data(:,index);
        else
            signal.(matlab.lang.makeValidName(tokens{i}{1}))=0;
        end
    end
    
    if torontoFile
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % identification of data for mass spectrometer related data sets
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
        signal.ts       =   signal.Time;        
        signal.Iv       =  -signal.Flow/1000;             % scaling Iv from l to m^3 (has other convention compared to Spiroware)      
        signal.He       =   signal.He/100;                % scaling % to 1
        signal.SF6      =   signal.SF6/100;               % scaling % to 1
        signal.CO2      =   signal.CO2/100;               % scaling % to 1

        if isempty(signal.ts) || isempty(signal.Iv) || isempty(signal.He) || isempty(signal.SF6) || isempty(signal.CO2)
            error('Toronto File is missing a signal'); 
        end
        
        signal.Ivss     =   zeros(size(signal.ts));         % artificial Ivss signal (to avoid conditional statements later)        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optional smoothing of raw data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        signal.IvssSmooth   =   signal.Ivss;
        signal.IvssSmooth1  =   signal.Ivss;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Graphical output
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if parameters.Simulation.verb                          % optional plotting of the manoeuvre
            plotMolarFractions(figureNo, signal, file); 
        end         
    else
        errorFlag       =   0;
        signal.ts       =   signal.Time/1000;                 % scaling ms to s
        signal.Iv       =   signal.Flow/1000;                 % scaling Iv from l to m^3
        signal.MMss     =   signal.MMss/1000;                 % scaling MMss signal
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % input file type dependent data evaluation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        signal.CO2      =   signal.CO2/100; 
        signal.O2       =   signal.O2/100;                        % scaling % to 1
        
        if sum(signal.ts)==0 || sum(signal.Iv)==0 || sum(signal.MMss)==0 || sum(signal.CO2)==0 || sum(signal.O2)==0
            if MissPlexy==1 && sum(signal.CO2)==0 
                errorFlag=0;
            else
                errorFlag=1;    
            end
        end
        
        
        if sbwPercent==1    %single breath washout
            if errorFlag==1 || sum(signal.MMms)==0 
                error('File is missing a signal (sbwPercent)');
            end
            
            signal.MMms     =   signal.MMms/1000;                       % scaling g/mol to kg/mol
            signal.Ivss     =  -ones(size(signal.ts))*sampleFlow;       % artificial sample flow signal
        elseif dtgmbw   
            if aFile
                if errorFlag==1 || sum(signal.MMms)==0
                    error('File is missing a signal (dtgmbw)');
                end
                %error('Doubletracer Gas Multiple Breath Washout not implemented for A Files');
                signal.MMms     =   signal.MMms/1000;                 % scaling g/mol to kg/mol
                signal.Ivss     =  -ones(size(signal.ts))*sampleFlow; 	% artificial sample flow signal
                
                signal.SF6      =   zeros(size(signal.MMms));
            else             
                if errorFlag==1 || sum(signal.N2)==0 || sum(signal.SF6)==0 || sum(signal.Ar)==0 || sum(signal.MMms)==0
                    error('File is missing a signal (dtgmbw)');
                end

                signal.N2       =   signal.N2/100;                    % scaling % to 1
                signal.SF6      =   signal.SF6/100;                   % scaling % to 1
                signal.Ar       =   signal.Ar/100;                    % scaling % to 1
                signal.MMms     =   signal.MMms/1000;                 % scaling g/mol to kg/mol
                signal.Ivss     =  -ones(size(signal.ts))*sampleFlow; 	% artificial sample flow signal

                signal.IvssSmooth   =   signal.Ivss;
                signal.IvssSmooth1  =   signal.Ivss;
            end
            
        elseif sbwMm==1   %Single breath washout      
            if errorFlag==1 || sum(signal.CO2)==0 
                error('File is missing a signal (sbwMm)');
            end
            
            signal.CO2      =   (signal.CO2-mmAir)/mm2Percent;          % evaluation of %-baesd CO2 signal
            signal.MMms     =   zeros(size(signal.ts));                 % dummy MMms signal
            signal.Ivss     =  -ones(size(signal.ts))*sampleFlow;       % artificial sample flow signal
            
        elseif heMbw==1     %Helium multiple breath washoud            
            if errorFlag==1  
                error('File is missing a signal (heMbw)');
            end
            
            signal.MMms     =   signal.MMms/1000;                       % scaling g/mol to kg/mol
            signal.Ivss     =  -ones(size(signal.ts))*sampleFlow;       % artificial sample flow signal
            
        elseif aFile==1 && spirowareType==0       
            if errorFlag==1 || sum(signal.MMms)==0 || sum(signal.SampleFlow)==0
                error('File is missing a signal (aFile and not a spiroware Type)');
            end
            
            signal.MMms     =   signal.MMms/1000;                     % scaling g/mol to kg/mol
            signal.Ivss     =   signal.SampleFlow/1000;                % scaling l/s to m^3/s
        elseif aFile==1
            if errorFlag==1 || sum(signal.MMms)==0 || sum(signal.SampleFlow)==0
                error('File is missing a signal (aFile)');
            end
            if sum(signal.MMms) == 0
                signal.MMms	=   ones(size(signal.ts))*(28.85/1000);     % artificial MMms signal
            else
                signal.MMms	=   signal.MMms/1000;             % scaling g/mol to kg/mol
            end
            if sum(signal.SampleFlow) == 0
                signal.Ivss	=  -ones(size(signal.ts))*sampleFlow;       % artificial sample flow signal
            else
                signal.Ivss	=   signal.SampleFlow/1000;        % scaling l/s to m^3/s
            end
        elseif bFile==1
            signal.N2       =   signal.N2/100; 
            if sum(signal.SampleFlow) == 0
                signal.Ivss	=  -ones(size(signal.ts))*sampleFlow;       % artificial sample flow signal
            else
                signal.Ivss	=   signal.SampleFlow/1000;        % scaling l/s to m^3/s
            end
        end      
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Correcting for 100% O2 if MBW (that is, if >99% is found at all)
        % Algorithm: find sections with O2 level > 99% and use the median of the
        % resulting sample set to determine the tentative 100% correction factor.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if aFile == 1
            signal = correctO2Signal(signal);               
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Optional smoothing of raw data
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if aFile
%             signal = filterRawSignal(signal, nOrder, filterCutoff, dt); 
%             signal = filterDelayIV(signal,parameters);
%             signal = dataDelay(signal, parameters);
        else %Copy data to avoid problems in later code
            signal.IvssSmooth   =   signal.Ivss;
            signal.IvssSmooth1  =   signal.Ivss;
        end
     
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Graphical output
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if parameters.Simulation.verb                           % optional plotting of the manoeuvre            
            plotOptionalDataFiltering(file, figureNo, signal);            
        end        
    end
    
end









%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% read file and split into header and data section
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Header,Data] = readFile(file)
    
    fid = fopen(file, 'r');
    if fid<1
        error('file %s not found\n',file);
    end
    
    Header = fgetl(fid);
    fclose(fid);
    
    %The csvread does not seem to work in the 2014b version of Matlab
    %Data   = csvread(file, 1, 0);  
    Data   = dlmread(file, '\t', 1,0);
end

function data = reshapeData(aFile, token, bFile, data, file, parameters)

    dtgmbw          =   parameters.Simulation.dtgmbwAnalysis;   %General dtg measurement
    dtgmbwBaby      =   parameters.Simulation.dtgmbwBaby;       %4% DTG measurement

    if aFile==1 && strncmp('(SPW_V3.1',token(end),9)
        nt = length(token)-1;                           % header contains software version string (new A file)
    elseif bFile==1
        nt = length(token)-1;
    elseif aFile==1 && dtgmbw && strncmp('(SPW_V3.2.1',token(end),11)
        nt = length(token)-1;  
    elseif aFile==1 && dtgmbw && strncmp('(SPW_V3.2.0',token(end),11)
        nt = length(token)-3;  
    elseif aFile==1 && dtgmbw
       nt = length(token)-3;     
    elseif aFile==1 && strncmp('(SPW_V3.2',token(end),9)
        nt = length(token)-3;      
    else
        nt = length(token);
    end
    ndata           =   round(length(data)/nt);         % number of data rows
        
    if nt*ndata ~= length(data)
        fprintf('nt=%d,ndata=%d,length(data)=%d',nt,ndata,length(data))
        error('readSignalSpirowareMassSpec: data inconsistent in file <%s>',file)
    end
    data            =   reshape(data,nt,ndata)';        % reshape accordingly
end




        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Correcting for 100% O2 if MBW (that is, if >99% is found at all)
% Algorithm: find sections with O2 level > 99% and use the median of the
% resulting sample set to determine the tentative 100% correction factor.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function signal = correctO2Signal(signal)
    meanO2=1;
    %     stdO2=.01;    % data from Zürich was not corrected with a to low threshold value
    stdO2=.05;
    discriminantO2 = signal.O2(signal.O2>(meanO2-stdO2));
    if length(discriminantO2)>10
        meanO2=median(discriminantO2);
        if meanO2<1
            meanO2=1.00075*meanO2; % empirical factor to correct for underestimation of O2 level
        elseif meanO2>1
            meanO2=0.99925*meanO2;
        end
        fprintf('O2 corrected by factor %g to correct 100%% O2 level\n',round(1/meanO2*1000)/1000);
    else
        fprintf('O2 was not corrected because the maximum measured O2 level was: %g \n',max(signal.O2)*100);
    end
end


       
function plotMolarFractions(figureNo, signal, file) 
    figure(200+figureNo);
    subplot(2,1,1);
    Vol  =   cumtrapz(signal.ts,signal.Iv); % integration of volume (plain)
    plot(signal.ts,signal.delayLine,signal.ts,signal.Iv*1000,signal.ts,Vol*1000);
    title(['file: "' file '"'],'Interpreter','none')
    xlabel('t / s')
    ylabel('delay / -, Iv / l/s, V / l')
    legend('delay', 'Iv mean', 'V_p_l_a_i_n')

    subplot(2,1,2);
    plot(signal.ts,signal.He*100,'b',signal.ts,signal.SF6*100,'r',signal.ts,signal.CO2*100,'g');
    title('Molar fractions')
    xlabel('t / s')
    ylabel('x / %')
    legendStrings={'He', 'SF6', 'CO2'};
    legend(legendStrings)
end

function plotOptionalDataFiltering(file, figureNo, signal) 
    [filePath,fileName,fileExtention] = fileparts(file);
    if fileName(1) == 'B'
        addFigure=10000;
    else
        addFigure=0;
    end

    figure(addFigure+100+figureNo)
    subplot(2,1,1);
    Vol  =   cumtrapz(signal.ts,signal.Iv); % integration of volume (plain)
    plot(signal.ts,signal.delayLine,signal.ts,signal.Iv*1000,signal.ts,Vol*1000);
    title('delay line and flow!')
    xlabel('t / s')
    ylabel('delay / -, Iv / l/s, V / l')
    legend('delay','Iv raw data','Vol raw data')

    subplot(2,1,2);
    plot(signal.ts,signal.delayTemp);
    title('delay temp!')
    xlabel('t / s')
    ylabel('temp / °C')
    legend('tempDelay')

    figure(200+figureNo);
    subplot(3,1,1);

    Vol  =   cumtrapz(signal.ts,signal.Iv); % integration of volume (plain)
    plot(signal.ts,signal.Iv*1000,signal.ts,Vol*1000);
    title(['file: "' file '"'],'Interpreter','none')
    xlabel('t / s')
    ylabel('Iv, V / l/s, l')
    legend('Iv mean', 'V_p_l_a_i_n')

    subplot(3,1,2);
    plot(signal.ts,signal.CO2*100,'b',signal.ts,signal.O2*100,'g');
    title('Molar fractions')
    xlabel('t / s')
    ylabel('x / %')
    legendStrings={'CO2', 'O2'};
    legend(legendStrings)

    subplot(3,1,3);
    plot(signal.ts,signal.MMss*1000,'b',signal.ts,signal.MMms*1000,'g');
    title('Molar masses')
    xlabel('t / s')
    ylabel('MMms / g/mol')
    legendStrings={'MMss', 'MMms'};
    legend(legendStrings) 
end
