%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% MATLAB function to transform file name in anonymous form
%
% Copy right: NM Numerical Modelling GmbH
% This model must not be distributed without explicit consent by NM GmbH
%
% Version 2.13, 10. Mai 2012
% Markus Roos, LuFu, Inselspital Bern
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [nameAnonymous,hash]=nameAnonymizer(name,hash)

    positions = strfind(name,'-');
    if length(positions) ~= 4
        warning('file name <%s> cannot be anonymized due to wrong structure\n',name);
    end
    
    dateExperiment=name([positions(1)+1:positions(2)-1]);
    timeExperiment=name([positions(2)+1:positions(3)-1]);
    person=name([positions(3)+1:positions(4)-1]);
    hashMatch=[dateExperiment,'-',person];
    
    i=1;
    [n,m]=size(hash);
    while i<=n
        if strcmp(hashMatch,hash{i,1})
            break;
        end
        i=i+1;
    end
    j=1;
    if i>n
        hash{i,1}=hashMatch;
        hash{i,2}={timeExperiment};
    else
        times=hash{i,2};
        while j<=length(times)
            if strcmp(timeExperiment,times{j})
                break;
            end
            j=j+1;
        end
        if j>length(times)
            hash{i,2}={times{:},timeExperiment};
        end
    end
    nameAnonymous=[name(1:positions(1)),sprintf('%d.%d',i,j),name([positions(1):positions(3)-1]),name(positions(4):length(name))];
end

