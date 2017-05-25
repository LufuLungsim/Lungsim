%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optional evaluation of Capno coefficients
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 12. Dezember 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [capno] = evaluateCapnoIndices(table,parameters)

    capnoCheck=parameters.Simulation.CapnoCheck;

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of capno indices (statictics)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    volExp   	= table.volExp(:,1);
    [~,index]   = sort(volExp,1,'ascend');
    nBreaths    = length(volExp);
    firstBreath = ceil( nBreaths*0.025);
    lastBreath  = floor(nBreaths*0.975);
    validIndices= index(firstBreath:lastBreath);
    okIndices   = find(cellfun(@(x) x.isValid==1, table.CO2.capno));%table.CO2.capno.isValid==1);
                        
    capno.validIndices=intersect(okIndices,validIndices);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of mean and std values
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    capno=evaluateMeanStd(capno,table.CO2.capno);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interactive data selection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if capnoCheck
        capnoArray = table.CO2.capno(:);
        [hfig,subPlotHandle]=plotCapnoIndices(table,capno,[0 1 1],901);
        while waitforbuttonpress()==0
            [xMin,xMax,yMin,yMax,axesClicked]=selectRubberband(gca);

            switch axesClicked
                case subPlotHandle(1)
                    yValues=cellfun(@(x) x.fowlerDead,capnoArray)*1e6;
                case subPlotHandle(2)
                    yValues=cellfun(@(x) x.coeffII,capnoArray)/1e3;
                case subPlotHandle(3)
                    yValues=cellfun(@(x) x.coeffIII,capnoArray)/1e3;
                case subPlotHandle(4)
                    yValues=cellfun(@(x) x.fCO2E,capnoArray)*1e2;
                case subPlotHandle(5)
                    yValues=cellfun(@(x) x.CO2et,capnoArray)*1e2;
                otherwise
                    yValues=[]; % not in correct window
            end
            if ~isempty(yValues)
                indexSelected=find(volExp>xMin/1e6 & volExp<xMax/1e6 & yValues>yMin & yValues<yMax);
                indexReselected=setdiff(indexSelected,capno.validIndices);
                capno.validIndices=setdiff(capno.validIndices,indexSelected);
                capno.validIndices=union(capno.validIndices,indexReselected);
                capno=evaluateMeanStd(capno,table.CO2.capno);
                [hfig,subPlotHandle]=plotCapnoIndices(table,capno,[0 1 1],901);
            end
        end
        close(hfig);
    end
end


function [xMin,xMax,yMin,yMax,axesClicked] = selectRubberband(handle)
% return rubberband coordinates

    point1 = get(handle, 'CurrentPoint'); % Button down detected
    rbbox(); % Return figure units
    point2 = get(handle, 'CurrentPoint'); % Button up detected
    axesClicked=gca;
    
    % Get rbbox limits
    xMin = min(point1(1, 1), point2(1, 1));
    xMax = max(point1(1, 1), point2(1, 1));
    yMin = min(point1(1, 2), point2(1, 2));
    yMax = max(point1(1, 2), point2(1, 2));

end

function [capnoIndices] = evaluateMeanStd(capnoIndices,CO2capno)
% return mean and std deviations of capno indices

    if ~isempty(capnoIndices.validIndices)
        capno= CO2capno(capnoIndices.validIndices);
        capnoIndices.meanfCO2E    	= 	mean(cellfun(@(x) x.fCO2E, capno));   
        capnoIndices.meanCO2et      =   mean(cellfun(@(x) x.CO2et, capno));
        capnoIndices.meanCoeffII    =   mean(cellfun(@(x) x.coeffII, capno));
        capnoIndices.meanCoeffIII   =   mean(cellfun(@(x) x.coeffIII, capno));
        capnoIndices.meanFowlerDead	=   mean(cellfun(@(x) x.fowlerDead, capno));
        
        capnoIndices.stdfCO2E       = 	std(cellfun(@(x) x.fCO2E, capno));
        capnoIndices.stdCO2et       =   std(cellfun(@(x) x.CO2et, capno));
        capnoIndices.stdCoeffII     =   std(cellfun(@(x) x.coeffII, capno));
        capnoIndices.stdCoeffIII    =   std(cellfun(@(x) x.coeffIII, capno));
        capnoIndices.stdFowlerDead	=   std(cellfun(@(x) x.fowlerDead, capno));
    else
        capnoIndices.meanfCO2E    	= 	0;
        capnoIndices.meanCO2et      = 	0;
        capnoIndices.meanCoeffII    = 	0;
        capnoIndices.meanCoeffIII   = 	0;
        capnoIndices.meanFowlerDead	= 	0;
        
        capnoIndices.stdfCO2E       = 	0;
        capnoIndices.stdCO2et       = 	0;
        capnoIndices.stdCoeffII     = 	0;
        capnoIndices.stdCoeffIII    = 	0;
        capnoIndices.stdFowlerDead	= 	0;
    end
    
end
