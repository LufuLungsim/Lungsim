%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optional evaluation of Scond/Sacin for given gas species
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.0, 08. Jan. 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [gas] = evaluateScondSacin(choice,lBound,uBound,gas,table,parameters)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Evaluation of Scond/Sacin and data to plot the corresponding figure
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exclLimit       =   parameters.Simulation.exclusionLimit;
    onlySmall       =   parameters.Simulation.onlySmall;
    ScondSacinCheck =   parameters.Simulation.ScondSacinCheck;
    ScondSacinNorm  =   parameters.Simulation.ScondSacinNorm;
    
    if ~isfield(gas.General,'scondSacin')
        gas.General.scondSacin={};
    end
    
    turnovers       = gas.General.to;
    volExp          = table.volExp(:,1);
    volExpMean      = mean(table.volExp);
    
    if ScondSacinNorm
        slopeType   = gas.slope(:,2)/1000;
    else
        slopeType   = gas.slope(:,2).*volExp;
    end
    
    % samples selection for Scond
    alpha = 1 - volExp/volExpMean; % alpha = abs(alpha); % alpha<0 means volExp>volExpMean which is admissible!
    if ~onlySmall
        alpha = abs(alpha);
    end
    indices.Used     = find(lBound < turnovers & uBound > turnovers & exclLimit > alpha);
    indices.Excluded = find(lBound < turnovers & uBound > turnovers & exclLimit <= alpha);
    indices.Possible = union(indices.Used,indices.Excluded);
    
    gas.General.scondSacin=setScondSacin(gas.General.scondSacin,choice,indices,turnovers,slopeType,lBound,uBound);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Interactive data selection
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    hDetail=0;
    if ScondSacinCheck
        [hfig]=plotScondSacin(choice,table,gas,parameters,[0 1 1],9011);
        while waitforbuttonpress()==0
            [xMin,xMax,yMin,yMax]=selectRubberband(gca);
            if xMin==xMax && yMin==yMax
                hfigPosition=get(hfig,'Position');
                hDetail=showDetails(xMin,yMin,turnovers,slopeType,indices,gas,1911,hfigPosition);
                figure(hfig);
                continue;
            end
            if ~isempty(turnovers)
                indicesSelected=find(turnovers>xMin & turnovers<xMax & slopeType>yMin & slopeType<yMax);
                indicesReselected=intersect(indices.Possible,setdiff(indicesSelected,indices.Used));
                indicesSelected=intersect(indicesSelected,indices.Possible);
                indicesSelected=setdiff(indices.Used,indicesSelected);
                indices.Used=union(indicesSelected,indicesReselected);
                gas.General.scondSacin=setScondSacin(gas.General.scondSacin,choice,indices,turnovers,slopeType,lBound,uBound);
                [hfig]=plotScondSacin(choice,table,gas,parameters,[0 1 1],9011);
            end
        end
        close(hfig);
        if hDetail~=0
            close(hDetail);
        end
    end   
end



function [hfig] = showDetails(xMin,yMin,turnovers,slopeType,indices,gas,figureNo,position)
% show breath details for graphical inspection
    limitDistance=0.06;

    minDistance=1e99;
    for j=1:length(indices.Possible)
        i=indices.Possible(j);
        distance=norm([xMin-turnovers(i);yMin-slopeType(i)]);
        if distance<minDistance
            minDistance=distance;
            iSelected=i;
        end
    end
    if minDistance<limitDistance
        volume=gas.breaths(iSelected).Volume;
        flow=gas.breaths(iSelected).Iv;
        species=gas.breaths(iSelected).Species;
        
        valueMax=max(volume);
        valueMin=min(volume);
        
        nSpecies=size(species,2);
    
        hfig=figure(figureNo);
        hfigPosition=[position(1)+position(3)+20,position(2),position(3),position(4)];
        set(hfig,'Position',hfigPosition);
        set(hfig,'Name',strcat('ScondSacin visualization of breath no=', num2str(iSelected)));
        set(hfig,'NumberTitle','off');
        subplot(2,1,2)
        plot(volume,flow,'b');
        title('volume flow plot');
        xlabel('Vol / m^3');
        ylabel('Iv / m^3/');
        axis([valueMin,valueMax,-inf,inf]);
        subplot(2,1,1)
        hold off
        plot(volume,species(:,1),'b');
        if nSpecies>1
            hold on
            plot(volume,species(:,2),'g');
        end
        title('Species signal plot');
        xlabel('Vol / m^3');
        ylabel('signal');
        axis([valueMin,valueMax,min(min(species)),max(max(species))]);
    else
        hfig=0;
        fprintf('nothing selected\n')
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

%choice 1 = scond/sacin
%choice 2 = scond*/sacin*
function [ScondSacin] = setScondSacin(ScondSacin,choice,indices,turnovers,slopeType,lBound,uBound)
% return mean and std deviations of capno indices
    if ~isempty(indices.Used) && choice == 2 && indices.Used(1)==1    %If choice == Scond* and the first indices is the first breath
        indices.Used(1)=[]; 
    end
    
    if ~isempty(indices.Used)
        indices.Excluded=setdiff(indices.Possible,indices.Used);

        % Scond calculation by means of linear polynomial least square fit      
        xValues  = turnovers(indices.Used);
        yValues  = slopeType(indices.Used);
        
        
        fitScond = polyfit(xValues,yValues,1);
        Scond = fitScond(1);

        residSlopes = yValues - polyval(fitScond,xValues);
        SSresid = sum(residSlopes.^2);
        SStotal = (length(yValues)-1) * var(yValues);
        r2      = 1-SSresid/SStotal;

        % Sacin calculation
        sacin = Scond*turnovers;
        sacin = slopeType-sacin;
        Sacin = sacin(1);
   
        ScondSacin.Scond{choice} = Scond;
        ScondSacin.Sacin{choice} = Sacin;
        ScondSacin.InclusionIndices{choice} = indices.Used;
        ScondSacin.ExclusionIndices{choice} = indices.Excluded;
        ScondSacin.Fit{choice} = fitScond;
        ScondSacin.R2{choice} = r2;
        ScondSacin.lBound{choice} = lBound;
        ScondSacin.uBound{choice} = uBound;
        ScondSacin.xValues{choice} = xValues;
        ScondSacin.yValues{choice} = yValues;
    else
        ScondSacin.Scond{choice} = 0;
        ScondSacin.Sacin{choice} = 0;
        ScondSacin.InclusionIndices{choice} = [];
        ScondSacin.ExclusionIndices{choice} = [];
        ScondSacin.Fit{choice} = 0;
        ScondSacin.R2{choice} = 0;
        ScondSacin.lBound{choice} = 0;
        ScondSacin.uBound{choice} = 0;
        ScondSacin.xValues{choice} = [];
        ScondSacin.yValues{choice} = [];
    end
    
end

