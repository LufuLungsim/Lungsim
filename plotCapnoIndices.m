%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Optional plotting ScondSacin data (with given color scheme)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.0, 12. Dezember 2014
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hfig,subPlotHandle] = plotCapnoIndices(table,capnoMean,color,figureNumber)

    if isempty(capnoMean.validIndices)
        warning('no capno index data left')
        hfig=figure(figureNumber);
        subPlotHandle=[subplot(5,1,1),subplot(5,1,2),subplot(5,1,3),subplot(5,1,4),subplot(5,1,5)];
        return;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plotting capno indices
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    indices    = capnoMean.validIndices;
    capnoValid = table.CO2.capno(indices);
    capnoAll   = table.CO2.capno;
    vol        = table.volExp(indices,1);
    
    hfig=figure(figureNumber);
    if figureNumber>900
        set(hfig,'Name','Capno indices: user assisted data selection');
    else
        set(hfig,'Name','Capno indices: graphical representation');
    end
    figureSize=get(hfig,'Position');
    height=figureSize(4);
    width=figureSize(3);
    if width>height
        figureSize(4)=1.9*height;
        figureSize(2)=figureSize(2)-0.9*height;
        set(hfig,'Position',figureSize);
    end
    fowlerDeadAll  =cellfun(@(x) x.fowlerDead,capnoAll  );
    fowlerDeadValid=cellfun(@(x) x.fowlerDead,capnoValid);
    subPlotHandle(1)=subplot(5,1,1);
    hold off
    plot(table.volExp(:,1)*1e6,fowlerDeadAll*1e6,  'ok','MarkerEdgeColor',color,'MarkerFaceColor',[1 1 1],'MarkerSize',5)
    hold on
    plot(vol*1e6,              fowlerDeadValid*1e6,'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,  'MarkerSize',5)
    plot([min(vol),max(vol)]*1e6,[capnoMean.meanFowlerDead,capnoMean.meanFowlerDead]*1e6,'k-.')
    axis([-inf,inf,0,1.5*max(fowlerDeadAll*1e6)]);
    xlabel('V_{exp} / ml'); ylabel('V_{Fowler} / ml')
    if figureNumber>900
        title({'<mouse down/up> shows rubber band to select points, <enter> accepts choice','','distribution of CO2-based Fowler dead space'})
    else
        title('distribution of CO2-based Fowler dead space')
    end
    hold off
    
    coeffIIAll  =cellfun(@(x) x.coeffII,capnoAll  );
    coeffIIValid=cellfun(@(x) x.coeffII,capnoValid);
    subPlotHandle(2)=subplot(5,1,2);
    hold off
    plot(table.volExp(:,1)*1e6,coeffIIAll/1e3,  'ok','MarkerEdgeColor',color,'MarkerFaceColor',[1 1 1],'MarkerSize',5)
    hold on
    plot(vol*1e6,              coeffIIValid/1e3,'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,  'MarkerSize',5)
    plot([min(vol),max(vol)]*1e6,[capnoMean.meanCoeffII,capnoMean.meanCoeffII]/1e3,'k-.')
    axis([-inf,inf,-inf,inf]);
    axis 'auto y'
    xlabel('V_{exp} / ml'); ylabel('S_{nII} / 1/L')
    title('CO2-based normalized phase II slopes (SnII)')
    hold off
    
    coeffIIIAll  =cellfun(@(x) x.coeffIII,capnoAll  );
    coeffIIIValid=cellfun(@(x) x.coeffIII,capnoValid);
    subPlotHandle(3)=subplot(5,1,3);
    hold off
    plot(table.volExp(:,1)*1e6,coeffIIIAll/1e3,  'ok','MarkerEdgeColor',color,'MarkerFaceColor',[1 1 1],'MarkerSize',5)
    hold on
    plot(vol*1e6,              coeffIIIValid/1e3,'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,  'MarkerSize',5)
    plot([min(vol),max(vol)]*1e6,[capnoMean.meanCoeffIII,capnoMean.meanCoeffIII]/1e3,'k-.')
    axis([-inf,inf,-inf,inf]);
    axis 'auto y'
    xlabel('V_{exp} / ml'); ylabel('S_{nIII} / 1/L')
    title('CO2-based normalized phase III slopes (SnIII)')
    hold off
    
    fCO2EAll  =cellfun(@(x) x.fCO2E,capnoAll  );
    fCO2EValid=cellfun(@(x) x.fCO2E,capnoValid);
    subPlotHandle(4)=subplot(5,1,4);
    hold off
    plot(table.volExp(:,1)*1e6,fCO2EAll*1e2,  'ok','MarkerEdgeColor',color,'MarkerFaceColor',[1 1 1],'MarkerSize',5)
    hold on
    plot(vol*1e6,              fCO2EValid*1e2,'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,  'MarkerSize',5)
    plot([min(vol),max(vol)]*1e6,[capnoMean.meanfCO2E,capnoMean.meanfCO2E]*1e2,'k-.')
    axis([-inf,inf,0,8]);
    xlabel('V_{exp} / ml'); ylabel('F_{CO2,exp} / ml')
    title('distribution of expired CO2-volume')
    hold off
    
    CO2etAll  =cellfun(@(x) x.CO2et,capnoAll  );
    CO2etValid=cellfun(@(x) x.CO2et,capnoValid);
    subPlotHandle(5)=subplot(5,1,5);
    hold off
    plot(table.volExp(:,1)*1e6,CO2etAll*1e2,  'ok','MarkerEdgeColor',color,'MarkerFaceColor',[1 1 1],'MarkerSize',5)
    hold on
    plot(vol*1e6,              CO2etValid*1e2,'ok','MarkerEdgeColor',color,'MarkerFaceColor',color,  'MarkerSize',5)
    plot([min(vol),max(vol)]*1e6,[capnoMean.meanCO2et,capnoMean.meanCO2et]*1e2,'k-.')
    axis([-inf,inf,0,8]);
    xlabel('V_{exp} / ml'); ylabel('CO2et / %')
    title('distribution of end expiratory CO2 concentration')
    hold off
end
