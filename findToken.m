%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Finding entry of token in a list (cell array)
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 1.7, 21. Mai 2011
% Markus Roos, NM GmbH
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input: Token, Header containing Tokens, Search direction
function  [n] = findToken(t,list,direction)
    
    for j=1:length(t)
        nList = length(list);
        if direction==1
            for i=1:nList
                if( isequal(t{j},list{i} ))
                    n=i;
                    return
                end
            end
        else
            for i=nList:-1:1
                if( isequal(t{j},list{i}) )
                    n=i;
                    return
                end
            end
        end
    end 
    n=0;
end