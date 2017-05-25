%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optional plotting ScondSacin data (with given color scheme)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.0, 08. Jan. 2015
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hfig] = plotScondSacin(choice,table,gas,parameters,color,figureNumber)

    ScondSacin      =   parameters.Simulation.ScondSacin;
    ScondSacinNorm  =   parameters.Simulation.ScondSacinNorm;
    species         =   gas.name;
    
    if choice==2 
        species = strcat('Star',species);
    end

    if ScondSacin == 1
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Plotting Scond/Sacin
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        lBound  = gas.General.scondSacin.lBound{choice};
        uBound  = gas.General.scondSacin.uBound{choice};
        Scond   = gas.General.scondSacin.Scond{choice};
        Sacin   = gas.General.scondSacin.Sacin{choice};
        r2      = gas.General.scondSacin.R2{choice};
        
        r2Annotation    = ['R^2 = ' num2str(r2,2)];
        
        ind = gas.General.scondSacin.InclusionIndices{choice};
        hasIncludes = sum(ind);
        indExcl = gas.General.scondSacin.ExclusionIndices{choice};
        hasExcludes = sum(indExcl);
        
        % Scond linear fit
        toValues = (lBound:0.02:uBound);
        fitSlopes = polyval(gas.General.scondSacin.Fit{choice},toValues);
        
        slopes = gas.slope(:,2);
        if ScondSacinNorm
            slopesType = slopes/1000;
        else
            slopesType = slopes.*table.volExp(:);
        end
        
        % Plots
        hfig=figure(figureNumber);
        if figureNumber>9000
            set(hfig,'Name','Scond and Sacin evaluation: user assisted data selection');
        else
            set(hfig,'Name','Scond and Sacin evaluation');
        end
     
        plot(gas.General.to(:),slopesType(:),...
            'ok',...
            'MarkerEdgeColor',color*0.7,...
            'MarkerFaceColor',color*0.7,...
            'MarkerSize',2)
        
        hold on
        
        if hasIncludes
            plot(gas.General.to(ind),slopesType(ind),...
                'ok',...
                'MarkerEdgeColor',color,...
                'MarkerFaceColor',color,...
                'MarkerSize',5)
        end
        
        if hasExcludes 
            plot(gas.General.to(indExcl),slopesType(indExcl),...
                'ok',...
                'MarkerEdgeColor',color,...
                'MarkerFaceColor',[1 1 1],...
                'MarkerSize',5)
        end
        
        plot(toValues,fitSlopes,...
            'Color',color,...
            'LineStyle','--',...
            'LineWidth',2)
        
        text(2,max(fitSlopes),r2Annotation,'Color',color);
        
        axis([0,round(max(gas.General.to))+1,min(0,min(slopesType)),round(5*(0.4+max(slopesType)))/5])
        
        xlabel('Turnover');
        
        if ScondSacinNorm
            ylabel('Sn_{III} [1/l]');
        else
            ylabel('Sn_{III}*V_{T,exp} [-]');
        end
        
        if ScondSacinNorm
            if figureNumber>9000
                title({'<mouse down/up> shows rubber band to select points, <enter> accepts choice','',[['Scond' species ' = '],num2str(Scond,3), ' [1/l/TO]; ',['Sacin'  species ' = '],num2str(Sacin,3), ' [1/l]']})
            else
                title([['Scond' species ' = '],num2str(Scond,3), ' [1/l/TO]; ',['Sacin'  species ' = '],num2str(Sacin,3), ' [1/l]'])
            end
        else
            if figureNumber>9000
                title({'<mouse down/up> shows rubber band to select points, <enter> accepts choice','',[['Scond' species ' = '],num2str(Scond,3), ' [1/TO]; ',['Sacin'  species ' = '],num2str(Sacin,3), ' [-]']})
            else
                title([['Scond' species ' = '],num2str(Scond,3), ' [1/TO]; ',['Sacin'  species ' = '],num2str(Sacin,3), ' [-]'])
            end
            
        end
        
        if ScondSacinNorm
            part = '';
        else
            part = '*V_{T,exp}';
        end
        
        index = 1;
        L{index} = ['SnIII' species part]; index=index+1;
        if hasIncludes
            L{index} = ['SnIII' species part ' included'];  index=index+1;
        end
        if hasExcludes
            L{index} = ['SnIII' species part ' excluded'];  index=index+1;
        end
        L{index} = ['fit to SnIII' species part];
        legend(L,'FontSize',10);
        hold off
    end
end