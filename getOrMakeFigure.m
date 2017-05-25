function h = getOrMakeFigure(figureName)
    
    if ~isempty(figureName)
        h = findobj('Name',figureName);
        if isempty(h)
            h = figure('Name',figureName,'NumberTitle','off');            
        end
    else
        h = figure('Name',figureName,'NumberTitle','off');
    end
    figure(h)
    
end
