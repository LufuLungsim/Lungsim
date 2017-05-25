%
% Determine the end expiratory Gas concentration through different methods
%
%
function N2et = getMMee(MMeeMethod, Gas, indicesHalfBreath, flow, startPercentage, endPercentage)

    volume = abs(cumsum(flow));
    [valMaxVol, indexMaxVolume] = max(volume);
    
    startInterval = find(volume(1:indexMaxVolume) <= valMaxVol*startPercentage, 1, 'last');
    endInterval   = find(volume(1:indexMaxVolume) <= valMaxVol*endPercentage,   1, 'last');
    
    switch MMeeMethod
        case 'MMeeR98'
            startInterval = find(volume(1:indexMaxVolume) <= valMaxVol*0.95, 1, 'last');
            endInterval   = find(volume(1:indexMaxVolume) <= valMaxVol*0.98, 1, 'last');
            N2et=median(Gas(indicesHalfBreath(startInterval:endInterval)));          
        case 'MMeeRMean' 
            N2et=mean(Gas(indicesHalfBreath(startInterval:endInterval)));
        case 'MMeeRMedian'
             N2et=median(Gas(indicesHalfBreath(startInterval:endInterval)));
        case 'MMeeRFit'
            GasFit=Gas(startInterval:endInterval);
            polynomFitInsp=polyfit((1:size(GasFit))',GasFit,1);
            N2et=polyval(polynomFitInsp,Gas(indicesHalfBreath(end)));
        otherwise
            msgbox('No MMee load method selected')
    end
    %Old method where the end tidal concentration was calculated based on
    %time and not volume
%     switch MMeeMethod
%         case 'MMeeR98'
%             N2et=median(Gas(indicesHalfBreath(round(0.95*length(indicesHalfBreath)):round(0.98*length(indicesHalfBreath)))));          
%         case 'MMeeRMean'
%             N2et=mean(Gas(indicesHalfBreath(max(1, round(startInterval*length(indicesHalfBreath))):round(endInterval*length(indicesHalfBreath)))));
%         case 'MMeeRMedian'
%             N2et=median(Gas(indicesHalfBreath(max(1, round(startInterval*length(indicesHalfBreath))):round(endInterval*length(indicesHalfBreath)))));
%         case 'MMeeRFit'
%             N2fit=Gas(indicesHalfBreath(max(1, round(startInterval*length(indicesHalfBreath))):round(endInterval*length(indicesHalfBreath))));
%             polynomFitInsp=polyfit((1:size(N2fit))',N2fit,1);
%             N2et=polyval(polynomFitInsp,Gas(indicesHalfBreath(end)));
%         otherwise
%             msgbox('No MMee load method selected')
%     end
end