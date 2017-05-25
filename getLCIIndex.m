%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Find index in list
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [index]=getLCIIndex(indices,compatibility)
    
    index=0;
    notFound=1;
    i=1;
    if compatibility    % return the very first breath with ee below the critical value
        if ~isempty(indices) 
            index=indices(1);
        end
        return;
    end
    while (i<=length(indices)-2 && notFound)
        if ((indices(i+1)==indices(i)+1) && (indices(i+2)==indices(i)+2))
            index=indices(i);
            return;
        else
            i=i+1;
        end
    end

end