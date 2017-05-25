%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Reading breathing manoeuver data (A-files and older or wbreath export .txt-files
%
% Version 3.0, 25. Juli 2012
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [header,data,signal]=readSignalWbreath(file,figureNo,parameters)

    dataType        =   parameters.Simulation.SF6;      % SF6 or He
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % read file and reformat raw data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [header,data]   =   readFile(file);               	% reading data file
    token           =   splitToken(header);            	% header tags
    
    nt = length(token);
    ndata           =   round(length(data)/nt);         % number of data rows
    
    if nt*ndata ~= length(data)
        fprintf('nt=%d,ndata=%d,length(data)=%d\n',nt,ndata,length(data))
        error('readSignalWbreath:inconsistentData','data inconsistent in file <%s>\n',file)
    end
    data	=   reshape(data,nt,ndata)';        % reshape accordingly
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identification of data for mass spectrometer related data sets
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    iTid    =   findToken({'Time'},token,1);
    iVol    =   findToken({'Vol'},token,-1);
    iFlow   =   findToken({'Flow'},token,-1);
    iMM     =   findToken({'MM'},token,1);
    if dataType == 0    % He
        iMMss   =   findToken('MMss',token,1);
    end
    
    
    signal.ts       =   data(:,iTid)/1000;              % scaling ms to s
    signal.Iv       =   data(:,iFlow)/1000;             % scaling Iv from l to m^3 (has other convention compared to Spiroware)
    signal.Vol      =   data(:,iVol)/1000;              % scaling Iv from l to m^3 (has other convention compared to Spiroware)
    signal.MM       =   data(:,iMM)/1000;               % scaling g/mol to kg/mol
    if dataType == 0    % He
        signal.MMss =   data(:,iMMss)/1000;             % scaling g/mol to kg/mol
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % generate delay line
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % signal   =   delayIV(signal,parameters);
    % signal   =   delayTemp(signal,parameters);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Graphical output
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if parameters.Simulation.verb                          % optional plotting of the manoeuvre
        figure(1000+figureNo);
        subplot(2,1,1);
        plot(signal.ts,signal.Iv*1000,signal.ts,signal.Vol*1000);
        title(['file: "' file '"'],'Interpreter','none')
        xlabel('t / s')
        ylabel('Iv / l/s, V / l')
        legend('Iv', 'V_p_l_a_i_n')
        
        subplot(2,1,2);
        plot(signal.ts,signal.MM*1000,'b');
        if dataType == 0
            hold on
            plot(signal.ts,signal.MMss*1000,'r');
            hold off
        end
        title('Mean molar mass')
        xlabel('t / s')
        ylabel('MM / g/mol')
        if dataType == 0
            legend('MM','MMss')
        else
            legend('MM')
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
    Data   = fscanf(fid,'%f');
    fclose(fid);
    
end