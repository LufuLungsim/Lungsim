%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to breath table generation
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 3.2, 08. Jan. 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [table]=getBreathTableWbreath(breathIndex,signal,parameters)

    tidal           =   parameters.Simulation.tidalState || parameters.Simulation.tidalMeanState;
    SF6             =   parameters.Simulation.SF6;
    sideStream      =   parameters.Simulation.sideStream && (SF6 == 0);
    
    meanStateEE     =   parameters.Simulation.meanEE;
    lowerBoundEE    =   parameters.Simulation.lowerBoundEE;
    upperBoundEE    =   parameters.Simulation.upperBoundEE;
    
    meanStateEI     =   parameters.Simulation.meanEI;
    lowerBoundEI    =   parameters.Simulation.lowerBoundEI;
    upperBoundEI    =   parameters.Simulation.upperBoundEI;
    
    ts              =   signal.ts;
    vol             =   signal.Vol;
    Iv              =   signal.Iv;
    if sideStream
        MM          =   signal.MMss;
    else
        MM          =   signal.MM;
    end
    
    
    if vol(breathIndex(2)-1)>0
        table.startIndex = 1;                              % first half is inspiration
        table.nHalfBreaths = 2*floor((length(breathIndex)-1)/2);
    else
        table.startIndex = 2;                           	% first half is expiration
        table.nHalfBreaths = 2*floor((length(breathIndex)-2)/2);
    end
    
    for i=1:table.nHalfBreaths/2
        j=table.startIndex+2*(i-1);
        indicesInsp=(breathIndex(j  ):breathIndex(j+1)-1);
        indicesExp =(breathIndex(j+1):breathIndex(j+2)-1);
        table.VolInsp(i,1)= trapz(ts(indicesInsp),Iv(indicesInsp)); % evaluation of vol with given Iv (consistency with other Iv-based integrals)
        table.VolExp(i,1) =-trapz(ts(indicesExp ),Iv(indicesExp ));
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Method to determine N2ee by taking the median, median or extrapolated fit in user interval
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        chosenIndicesExp = indicesExp(-vol(indicesExp )>lowerBoundEE*table.VolExp(i,1) & -vol(indicesExp )<upperBoundEE*table.VolExp(i,1));
        if meanStateEE == 1
            table.MMee(i,1)=mean(MM(chosenIndicesExp));
        elseif meanStateEE ==2
            table.MMee(i,1)=median(MM(chosenIndicesExp));
        else    % mean state == 3
            polynomFitExp =polyfit(vol(chosenIndicesExp),MM(chosenIndicesExp),1);
            table.MMee(i,1)=polyval(polynomFitExp,vol(indicesExp (end)));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Method to determine N2ei by taking the median, median or extrapolated fit in user interval
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        chosenIndicesInsp = indicesInsp(vol(indicesInsp)>lowerBoundEI*table.VolInsp(i,1) & vol(indicesInsp)<upperBoundEI*table.VolInsp(i,1));
        if meanStateEI == 1
            table.MMei(i,1)=mean(MM(chosenIndicesInsp));
        elseif meanStateEE ==2
            table.MMei(i,1)=median(MM(chosenIndicesInsp));
        else    % mean state == 3
            polynomFitInsp=polyfit(vol(chosenIndicesInsp),MM(chosenIndicesInsp),1);
            table.MMei(i,1)=polyval(polynomFitInsp,vol(indicesInsp(end)));
        end
        
        table.VolMMInsp(i,1)            =   trapz(ts(indicesInsp),MM(indicesInsp).*Iv(indicesInsp));  % used for WO/WI detection
        table.VolMMExp(i,1)             =  -trapz(ts(indicesExp ),MM(indicesExp ).*Iv(indicesExp ));
        
        table.MMExpMean(i,1)            =	table.VolMMExp(i,1)/table.VolExp(i,1);
        table.MMInspMean(i,1)           =	table.VolMMInsp(i,1)/table.VolInsp(i,1);
        
        table.breathDuration(i,1)       =   ts(indicesExp(end))-ts(indicesInsp(1)-1);
        table.FlowInspMean(i,1)         =   mean(Iv(indicesInsp));
        table.FlowExpMean(i,1)          =  -mean(Iv(indicesExp ));
        
        if tidal
            tsInspStart                     =   ts(indicesInsp(1)-1);
            tsExpStart                      =   ts(indicesExp(1)-1);
            table.inspDuration(i,1)         =   ts(indicesInsp(end))-tsInspStart;
            table.expDuration(i,1)          =   ts(indicesExp(end))-tsExpStart;
            
            [table.inspPeakFlow(i,1),iMax]  =   max(Iv(indicesInsp));
            table.timeInspPeakFlow(i,1)     =   ts(indicesInsp(iMax))-tsInspStart;
            [table.expPeakFlow(i,1),iMax]	=   max(-Iv(indicesExp ));
            table.timeExpPeakFlow(i,1)      =   ts(indicesExp (iMax))-tsExpStart;
            
            table.VolTidal(i,1)             =   (table.VolInsp(i,1)+table.VolExp(i,1))/2;
            
            table.tPTEF2tERatio(i,1)        =	table.timeExpPeakFlow(i,1)/table.expDuration(i,1);
        end
        
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Additional evaluations
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    table.VolInspMean=mean(table.VolInsp);
    table.VolExpMean =mean(table.VolExp );
    table.VolInspStd =std(table.VolInsp,1)/table.VolInspMean;
    table.VolExpStd  =std(table.VolExp ,1)/table.VolExpMean;
    
    table.signal = signal;
    table.breathIndex = breathIndex;
    
end

