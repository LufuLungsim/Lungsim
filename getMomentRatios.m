% Calculates the moment ratios of a gas as described in 
% Moment ratio analysis of multiple breath nitrogen washout in infants with
% lung disease by Schibler et. al
% http://erj.ersjournals.com/content/15/6/1094.full.pdf
function [ m0, m1, m2 ] = getMomentRatios( gasConcentration, LCI )


%     m0 = sum(gasConcentration);
%     m1 = sum(gasConcentration.*LCI);
%     m2 = sum(gasConcentration.*LCI.^2);

    m0 = trapz(LCI, gasConcentration);
    m1 = trapz(LCI, gasConcentration.*LCI);
    m2 = trapz(LCI, gasConcentration.*LCI.^2);

end

