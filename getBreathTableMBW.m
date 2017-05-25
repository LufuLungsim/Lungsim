%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function for breath table generation
%
% The function is executed once per file and should transfer the
% information from the corrected signal to a table containing data about
% each breath. 
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ table ] = getBreathTableMBW( breathIndex,signal,parameters )
    
    %Initialisation
    n2mbw               =   parameters.Simulation.n2mbwAnalysis;                    % N2-MBW  (N2 Multiple Breath Washout)
    tidal               =   parameters.Simulation.tidalAnalysis;                    % tidal maneuver
    dtgmbw              =   parameters.Simulation.dtgmbwAnalysis;
    MissPlexy           =   parameters.Simulation.MissPlexy;
    
    torontoFile         =   parameters.Simulation.torontoFile;                      % based on mass spectrometer data
    
    interactiveNo       =   parameters.Simulation.interactiveNo;
    interactiveYes      =   parameters.Simulation.interactiveYes;
    interactiveDeviation=   parameters.Simulation.interactiveDeviation;
    interactive         =   interactiveNo+interactiveYes*2+interactiveDeviation*4;
    rmsCrit             =   parameters.Simulation.rmsCrit;
    bFile               =   parameters.Simulation.bFile;
    
    ScondSacin          =   parameters.Simulation.ScondSacin;
    ScondSacinCheck     =   parameters.Simulation.ScondSacinCheck;
    capnoIndices        =   parameters.Simulation.CapnoIndices;
    capnoCheck          =   parameters.Simulation.CapnoCheck;
    slopesLogFile       =   parameters.Simulation.slopesState;
    
    MMeeMethod          =   parameters.Simulation.MMeeMethod;
    MMeeTMin            =   parameters.Simulation.MMeeTMin;
    MMeeTMax            =   parameters.Simulation.MMeeTMax;
    
    
    DSpost              =   parameters.Device.volumePostcap;
    DSpre               =   parameters.Device.volumePrecap;
    btpsInsp            =   parameters.Calibration.factorBTPSinsp;
   
%     if MissPlexy==1;
%         btpsInsp        =   0;
%     end
%     if bFile == 1;
%        btpsInsp         =   0; 
%        DSpost           =   0;
%     end
    
    hasSF6              =   0;
    
    relativeMin=0.65;
    relativeMax=0.95;
    relativeCO2Min=0.05;
    relativeCO2Max=0.60;
    
    ts=signal.ts;
    vol=signal.VolBTEffBTPS;
    Iv=signal.IvEffBTPS;
   
    if torontoFile
        He= signal.He;
        SF6=signal.SF6;
        CO2=signal.CO2;
    else
        gases = {signal.N2Poi3, signal.O2Poi3, signal.CO2Poi3, signal.MMssPoi3};
        names = {'N2','O2','CO2','MMss'};
        MMss=signal.MMssPoi3;
        if dtgmbw && isfield(signal, 'SF6')
            gases(end+1) = {signal.SF6};
            names(end+1) = {'SF6'};
            hasSF6 = 1;
        end
    end

    table = getGasFromSignal(signal);
    
    %
    %Calculations
    %
    [table.startIndex, table.nHalfBreaths] = getHalfbreathSegmentation(vol, breathIndex);
    
    %% TODO
    % Make more abstract to treat all gases the same
    
    for i=1:table.nHalfBreaths/2
        j=table.startIndex+2*(i-1);
        indicesInsp             = breathIndex(j  ):breathIndex(j+1)-1;
        table.indicesInsp(i,1)  = indicesInsp(end);
        indicesExp              = breathIndex(j+1):breathIndex(j+2)-1;
        table.indicesExp(i,1)   = indicesExp(end);
        
        table.breathDuration(i,1)   = ts(indicesExp(end))-ts(indicesInsp(1)-1);
        table.volInsp(i,1)          = vol(indicesInsp(end));
        table.volExp(i,1)           =-vol(indicesExp(end) );
        table.flowInspMean(i,1)     = mean(Iv(indicesInsp));
        table.flowExpMean(i,1)      =-mean(Iv(indicesExp ));
        
        if capnoIndices
            if ~iscell(table.CO2.capno) && table.CO2.capno == 0
                table.CO2.capno = [];
            end
            if MissPlexy==1
                signal.CO2Poi3(indicesExp)  =   (1:length(indicesExp))*(0.05/length(indicesExp));
            end
            [table.CO2.capno{i,1},capnoCheck]=getCapno(relativeMin,relativeMax,relativeCO2Min,relativeCO2Max,-vol(indicesExp),signal.CO2Poi3(indicesExp),Iv(indicesExp),capnoCheck,rmsCrit);
        end
        
        if ~torontoFile
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % New method to treat all the gases equally
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            for j = 1:length(gases)
                table.(names{j}).et(i,1)            = getMMee(MMeeMethod, gases{j}, indicesExp,  Iv(indicesExp ), MMeeTMin, MMeeTMax);
                table.(names{j}).eti(i,1)           = getMMee(MMeeMethod, gases{j}, indicesInsp, Iv(indicesInsp), MMeeTMin, MMeeTMax); 
                table.(names{j}).exp(i,1)           =-trapz(ts(indicesExp ),gases{j}(indicesExp ).*Iv(indicesExp ));
                table.(names{j}).insp(i,1)          = trapz(ts(indicesInsp),gases{j}(indicesInsp).*Iv(indicesInsp));
                table.(names{j}).expRoof(i,1)       =-trapz(ts(indicesExp) ,(table.(names{j}).eti(i,1)-gases{j}(indicesExp)).*Iv(indicesExp ));
                table.(names{j}).inspRoof(i,1)      = trapz(ts(indicesInsp),(table.(names{j}).eti(i,1)-gases{j}(indicesInsp)).*Iv(indicesInsp ));                
                table.(names{j}).notReExp(i,1)      = DSpost*btpsInsp*(table.(names{j}).eti(i,1)-table.(names{j}).et(i,1)); 
                table.(names{j}).reInsp(i,1)        = DSpre*table.(names{j}).et(max(i-1,1))*btpsInsp;
            end
  
            
            table.MMss.slope(i,1:2)     = zeros(1,2);
            table.MMss.calcSlope(i,1:2) = zeros(1,2);
            table.MMss.diffSlope(i,1:2) = zeros(1,2);
            if slopesLogFile
                [table.CO2.slopePlain(i,1:2),interactive] = getSlope('lin','CO2',relativeMin,relativeMax,-vol(indicesExp),signal.CO2Poi3(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                [table.N2.slopePlain(i,1:2),interactive]  = getSlope('lin','N2',relativeMin,relativeMax,-vol(indicesExp),signal.N2Poi3(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                if ~hasSF6
                    [table.MMss.slope(i,1:2),interactive]     = getSlope('lin','MMss',relativeMin,relativeMax,-vol(indicesExp),MMss(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                    [table.MMss.calcSlope(i,1:2),interactive] = getSlope('lin','MMssCalc',relativeMin,relativeMax,-vol(indicesExp),signal.MMssCalc(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                    [table.MMss.diffSlope(i,1:2),interactive] = getSlope('lin','DiffMMss',relativeMin,relativeMax,-vol(indicesExp),signal.DiffMMss(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                end
            else
                table.CO2.slopePlain(i,1:2) = zeros(1,2);
                table.N2.slopePlain(i,1:2)  = zeros(1,2);
                table.SF6.slopePlain(i,1:2) = zeros(1,2);
            end
            table.CO2.slope(i,1:length(table.CO2.slopePlain(i,:))) = table.CO2.slopePlain(i,:);   % to simplify output routines
            
            if ScondSacin
                [table.N2.slope(i,1:2),interactive]=getSlope('normlin','N2',relativeMin,relativeMax,-vol(indicesExp),signal.N2Poi3(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                %%TODO Volumen / fluss etc in for schlaufe oben integrieren
                if ScondSacinCheck
                    table.N2.breaths(i,1).Volume=-vol(indicesExp);
                    table.N2.breaths(i,1).Iv=Iv(indicesExp);
                    table.N2.breaths(i,1).Species=signal.N2Poi3(indicesExp);
                end
                if hasSF6
                    [table.SF6.slope(i,1:2),interactive]=getSlope('normlin','SF6',relativeMin,relativeMax,-vol(indicesExp),signal.SF6(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                    if ScondSacinCheck
                        table.SF6.breaths(i,1).Volume=-vol(indicesExp);
                        table.SF6.breaths(i,1).Iv=Iv(indicesExp);
                        table.SF6.breaths(i,1).Species=signal.SF6(indicesExp);
                    end
                end
            else
                table.N2.slope(i,1:2)=zeros(1,2);
                table.SF6.slope(i,1:2)=zeros(1,2);
            end
            
        else
            table.He.et(i,1)     	= getMMee(MMeeMethod, He, indicesExp, Iv(indicesExp), MMeeTMin, MMeeTMax);
            table.He.reInsp(i,1) 	= DSpost*table.He.et(max(i-1,1))*btpsInsp;
            table.He.insp(i,1)   	= trapz(ts(indicesInsp),He(indicesInsp).*Iv(indicesInsp));
            table.He.exp(i,1)     	=-trapz(ts(indicesExp ),He(indicesExp ).*Iv(indicesExp ));
            table.He.meanInsp(i,1)  = mean(He(indicesInsp));
            
            table.SF6.et(i,1)     	= getMMee(MMeeMethod, signal.SF6, indicesExp, Iv(indicesExp), MMeeTMin, MMeeTMax);
            table.SF6.reInsp(i,1) 	= DSpost*table.SF6.et(max(i-1,1))*btpsInsp;
            table.SF6.insp(i,1)  	= trapz(ts(indicesInsp),signal.SF6(indicesInsp).*Iv(indicesInsp));
            table.SF6.exp(i,1)  	=-trapz(ts(indicesExp ),signal.SF6(indicesExp ).*Iv(indicesExp ));
            table.SF6.meanInsp(i,1)	= mean(signal.SF6(indicesInsp));
            
            if slopesLogFile
                [table.CO2.slopePlain(i,1:2),interactive] = getSlope('lin','CO2',relativeMin,relativeMax,-vol(indicesExp),CO2(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                [table.He.slopePlain(i,1:2),interactive]  = getSlope('lin','He', relativeMin,relativeMax,-vol(indicesExp),He(indicesExp), Iv(indicesExp),interactive,rmsCrit);
                [table.SF6.slopePlain(i,1:2),interactive] = getSlope('lin','SF6',relativeMin,relativeMax,-vol(indicesExp),signal.SF6(indicesExp),Iv(indicesExp),interactive,rmsCrit);
            else
                table.CO2.slopePlain(i,1:2) = zeros(1,2);
                table.He.slopePlain(i,1:2)  = zeros(1,2);
                table.SF6.slopePlain(i,1:2) = zeros(1,2);
            end
            
            if ScondSacin
                [table.CO2.slope(i,1:2),interactive]=getSlope('normlin','CO2',relativeMin,relativeMax,-vol(indicesExp),He(indicesExp), Iv(indicesExp),interactive,rmsCrit);
                [table.He.slope(i,1:2),interactive] =getSlope('normlin','He', relativeMin,relativeMax,-vol(indicesExp),He(indicesExp), Iv(indicesExp),interactive,rmsCrit);
                [table.SF6.slope(i,1:2),interactive]=getSlope('normlin','SF6',relativeMin,relativeMax,-vol(indicesExp),signal.SF6(indicesExp),Iv(indicesExp),interactive,rmsCrit);
                if ScondSacinCheck
                    table.He.breaths(i,1).Volume=-vol(indicesExp);
                    table.He.breaths(i,1).Iv=Iv(indicesExp);
                    table.He.breaths(i,1).Species=He(indicesExp);
                    table.SF6.breaths(i,1).Volume=-vol(indicesExp);
                    table.SF6.breaths(i,1).Iv=Iv(indicesExp);
                    table.SF6.breaths(i,1).Species=signal.SF6(indicesExp);
                end
            else
                table.CO2.slope(i,1:2)=zeros(1,2);
                table.He.slope(i,1:2)= zeros(1,2);
                table.SF6.slope(i,1:2)=zeros(1,2);
            end
        end
        
        %Dummy data to avoid errors
        table.InterceptionVolume(i,:)=0.0;          % not used for MBW
        table.DiffMMssMaxVolume(i,:)=0.0;           % not used for MBW
        table.EndTidalDiffMMss(i,:)=0.0;            % not used for MBW
        table.DiffMMssMax(i,:)=0.0;                 % not used for MBW
        
        if tidal || n2mbw || dtgmbw
            tsInspStart                     =   ts(indicesInsp(1)-1);
            tsExpStart                      =   ts(indicesExp(1)-1);
            table.inspDuration(i,1)         =   ts(indicesInsp(end))-tsInspStart;
            table.expDuration(i,1)          =   ts(indicesExp(end))-tsExpStart;

            [table.inspPeakFlow(i,1),iMax]  =   max(Iv(indicesInsp));
            table.timeInspPeakFlow(i,1)     =   ts(indicesInsp(iMax))-tsInspStart;
            [table.expPeakFlow(i,1),iMax]	=   max(-Iv(indicesExp ));
            table.timeExpPeakFlow(i,1)      =   ts(indicesExp (iMax))-tsExpStart;

            table.volTidal(i,1)             =   (table.volInsp(i,1)+table.volExp(i,1))/2;

            table.tPTEF2tERatio(i,1)        =	table.timeExpPeakFlow(i,1)/table.expDuration(i,1);
        end
    end

end

function [startIndex, nHalfBreaths] = getHalfbreathSegmentation(vol, breathIndex)
    maxVolume = max(vol(breathIndex(2):breathIndex(3)));
    minVolume = min(vol(breathIndex(2):breathIndex(3)));

    if abs(maxVolume)>abs(minVolume)
        startIndex = 2;                              % first half is inspiration
        nHalfBreaths = 2*floor((length(breathIndex)-2)/2);
    else
        startIndex = 3;                           	% first half is expiration
        nHalfBreaths = 2*floor((length(breathIndex)-3)/2);
    end
end

