%function to retrospectively find the point where the LCI was predicted
%with some certanity
function LCIStopCrit = getLCIFromPredictionStopCriteria(outputMeasurement, offsetGeneralInfo, lengthHeader, gas, predictionStructs, validData)
    
    %Output Measurement is a cell array with the fields as described in
    %generalInfo and header
    %generalInfo = {'Breaths used', 'Breaths excluded', 'ratio of Breaths used'}
    %header = {'LCI', 'Rel LCI Error', 'Breath 2_5', 'CEV', 'VN2Exp', 'FRC', 'R2 %/CEV', 'R2 %/#', 'R2 Vol/CEV', 'Estimated Error'};
    
    %Initialzie variables
    nBreaths            = gas.Standard.nBreaths(1);
    nBreaths            = nBreaths-sum(validData(1:nBreaths)==0);
    
    
    numOfMethods        = (size(outputMeasurement,2)-offsetGeneralInfo)/lengthHeader;
    nBreathsUsed        = cell2mat(outputMeasurement(3:end, 1));
    ratioOfBreaths      = cell2mat(outputMeasurement(3:end, 3));
    
    stopCrits           = {'stop3*1',  'stop3*0.75', 'stop3*0.5', 'stopN2005', 'stopTO6'};
    lengthStopCrits     = length(stopCrits);
    
    LCIStopCritVal      = zeros(numOfMethods*lengthStopCrits, 1);
    LCIStopCritRel      = zeros(numOfMethods*lengthStopCrits, 1);
    FRCStopCritVal      = zeros(numOfMethods*lengthStopCrits, 1);
    CEVStopCritVal      = zeros(numOfMethods*lengthStopCrits, 1);
    M0StopCritVal       = zeros(numOfMethods*lengthStopCrits, 1);
    M1StopCritVal       = zeros(numOfMethods*lengthStopCrits, 1);
    M2StopCritVal       = zeros(numOfMethods*lengthStopCrits, 1);
    LCIStopCritBreath   = zeros(numOfMethods*lengthStopCrits, 1);
    LCIStopCritRatio    = zeros(numOfMethods*lengthStopCrits, 1);
    
    for i = 1:numOfMethods
        pos = (i-1)*lengthHeader+offsetGeneralInfo;
        
        %Initialzize variables
        LCI             = cell2mat(outputMeasurement(3:end, pos+1));
        CEV             = cell2mat(outputMeasurement(3:end, pos+3)); 
        FRC             = cell2mat(outputMeasurement(3:end, pos+4)); 
        M0              = cell2mat(outputMeasurement(3:end, pos+5));        
        M1              = cell2mat(outputMeasurement(3:end, pos+6));       
        M2              = cell2mat(outputMeasurement(3:end, pos+7));
%         EndBreaths      = cell2mat(outputMeasurement(3:end, pos+3));        
%         CEV             = cell2mat(outputMeasurement(3:end, pos+4));
%         VN2Exp          = cell2mat(outputMeasurement(3:end, pos+5));
%         FRC             = cell2mat(outputMeasurement(3:end, pos+6));
        
        %Calculate relative difference to previous breath
        rdLCI           = (LCI(2:end)-LCI(1:end-1))./LCI(2:end);
%         rdEndBreaths    = (EndBreaths(2:end)-EndBreaths(1:end-1))./EndBreaths(2:end);
%         rdCEV           = (CEV(2:end)-CEV(1:end-1))./CEV(2:end);
%         rdVN2Exp        = (VN2Exp(2:end)-VN2Exp(1:end-1))./VN2Exp(2:end);
%         rdFRC           = (FRC(2:end)-FRC(1:end-1))./FRC(2:end);
        
        %Stop criteria where all the parameters dont change more than 5%
        %for 3 consecutive breaths
%         minRelDiff      = 0.05;
%         
%         indexLCI        = getPossibleEndpoints(rdLCI, minRelDiff);  
%         indexEndBreaths = getPossibleEndpoints(rdEndBreaths, minRelDiff);
%         indexCEV        = getPossibleEndpoints(rdCEV, minRelDiff);
%         indexVN2Exp     = getPossibleEndpoints(rdVN2Exp, minRelDiff);
%         indexFRC        = getPossibleEndpoints(rdFRC, minRelDiff);
%         
%         indexStop005Comb= find(indexLCI.*indexEndBreaths.*indexVN2Exp.*indexFRC.*indexCEV == 1, 1);
%         LCI005Comb      = LCI(indexStop005Comb);
%         breath005Comb   = nBreathsUsed(indexStop005Comb);
%         ratio005Comb    = ratioOfBreaths(indexStop005Comb);
        
        %Stop criteria where the LCI does not change more than 1% for 3
        %consecutive breaths
        indexLCI        = getPossibleEndpoints(rdLCI, 0.01);        
        indexStop0025   = find(indexLCI== 1, 1);
        
        if isempty(indexStop0025) || indexStop0025 > nBreaths
            indexStop0025 = nBreaths-5;
        end
        
        LCI0025         = LCI(indexStop0025);
        FRC0025         = FRC(indexStop0025);
        CEV0025         = CEV(indexStop0025);
        M0Stop0025      = M0(indexStop0025);
        M1Stop0025      = M1(indexStop0025);
        M2Stop0025      = M2(indexStop0025);
        breath0025      = nBreathsUsed(indexStop0025);
        ratio0025       = ratioOfBreaths(indexStop0025);
        
         %Stop criteria where the LCI does not change more than 0.75% for 3
        %consecutive breaths
        indexLCI        = getPossibleEndpoints(rdLCI, 0.0075);
        indexStop002    = find(indexLCI== 1, 1);
        
        if isempty(indexStop002) || indexStop002 > nBreaths
            indexStop002 = nBreaths-5;
        end
        
        LCI002         = LCI(indexStop002);
        FRC002         = FRC(indexStop002);
        CEV002         = CEV(indexStop002);
        M0Stop002      = M0(indexStop002);
        M1Stop002      = M1(indexStop002);
        M2Stop002      = M2(indexStop002);
        breath002      = nBreathsUsed(indexStop002);
        ratio002       = ratioOfBreaths(indexStop002);
        
        %Stop criteria where the LCI does not change more than 0.5% for 3
        %consecutive breaths
        indexLCI        = getPossibleEndpoints(rdLCI, 0.005);
        indexStop001    = find(indexLCI== 1, 1);
        
        if isempty(indexStop001) || indexStop001 > nBreaths
            indexStop001 = nBreaths-5;
        end
        
        LCI001         = LCI(indexStop001);
        FRC001         = FRC(indexStop001);
        CEV001         = CEV(indexStop001);
        M0Stop001      = M0(indexStop001);
        M1Stop001      = M1(indexStop001);
        M2Stop001      = M2(indexStop001);
        breath001      = nBreathsUsed(indexStop001);
        ratio001       = ratioOfBreaths(indexStop001);
        
       
        %Stop criteria where the lci does not change more than 5% for three
        %consecutive breaths and stays within a 5% band for 3 consecutive
        %breaths
%         indexLCI1       = getPossibleEndpoints(rdLCI, minRelDiff);  
%         indexLCI2       = getValAcordingToStopCriteria(rdLCI, minRelDiff);
%         
%         indexStop_3x3   = find(indexLCI1.*indexLCI2 == 1, 1);
%         LCI3x3Comb      = LCI(indexStop_3x3);
%         breath3x3Comb   = nBreathsUsed(indexStop_3x3);
%         ratio3x3Comb    = ratioOfBreaths(indexStop_3x3);
        
        
        
        %Get prediction for where N2 5% was reached
        breathN2005 = gas.Standard.nBreaths(4);
        indexN2005  = find(nBreathsUsed == breathN2005, 1, 'first');
        LCIN2005    = LCI(indexN2005);
        FRCN2005    = FRC(indexN2005 );
        CEVN2005    = CEV(indexN2005 );
        M0StopN2005 = M0(indexN2005);
        M1StopN2005 = M1(indexN2005);
        M2StopN2005 = M2(indexN2005);
        ratioN2005  = ratioOfBreaths(indexN2005);
        
        %Get prediction for where N2 9% was reached
        breathN2009 = gas.Standard.nBreaths(7);
        breathN009  = find(nBreathsUsed == breathN2009, 1, 'first');
        breathTO6   = find(gas.General.to>6, 1);
        indexTO6    = find(nBreathsUsed == breathTO6, 1, 'first');
        LCIN2009    = LCI(indexTO6);
        FRCN2009    = FRC(indexTO6);
        CEVN2009    = CEV(indexTO6);
        M0StopN2009 = M0(indexTO6);
        M1StopN2009 = M1(indexTO6);
        M2StopN2009 = M2(indexTO6);
        ratioN2009  = ratioOfBreaths(indexTO6);
        
       %Replace occurences where the stop criteria was not met with nan
                    [LCI0025,  breath0025, FRC0025, CEV0025, M0Stop0025,   M1Stop0025,   M2Stop0025,  ratio0025]  ...
       = nanIfEmpty( LCI0025,  breath0025, FRC0025, CEV0025, M0Stop0025,   M1Stop0025,   M2Stop0025,  ratio0025);
                    [LCI002,   breath002,  FRC002,  CEV002,  M0Stop002,    M1Stop002,    M2Stop002,   ratio002]   ...
       = nanIfEmpty( LCI002,   breath002,  FRC002,  CEV002,  M0Stop002,    M1Stop002,    M2Stop002,   ratio002);
                    [LCI001,   breath001,  FRC001,  CEV001,  M0Stop001,    M1Stop001,    M2Stop001,   ratio001]   ...
       = nanIfEmpty( LCI001,   breath001,  FRC001,  CEV001,  M0Stop001,    M1Stop001,    M2Stop001,   ratio001); 
                    [LCIN2005, breathN2005,FRCN2005,CEVN2005,M0StopN2005,  M1StopN2005,  M2StopN2005, ratioN2005] ...
       = nanIfEmpty( LCIN2005, breathN2005,FRCN2005,CEVN2005,M0StopN2005,  M1StopN2005,  M2StopN2005, ratioN2005);
                    [LCIN2009, breathN2009,FRCN2009,CEVN2009,M0StopN2009,  M1StopN2009,  M2StopN2009, ratioN2009] ...
       = nanIfEmpty( LCIN2009, breathN2009,FRCN2009,CEVN2009,M0StopN2009,  M1StopN2009,  M2StopN2009, ratioN2009);
        
        LCIStopCritVal(1+(i-1)*lengthStopCrits:i*lengthStopCrits)       = [LCI0025; LCI002; LCI001; LCIN2005; LCIN2009];
        LCIStopCritRel(1+(i-1)*lengthStopCrits:i*lengthStopCrits)       = ([LCI0025; LCI002; LCI001; LCIN2005; LCIN2009]-gas.Standard.lci(1))./gas.Standard.lci(1);
        CEVStopCritVal(1+(i-1)*lengthStopCrits:i*lengthStopCrits)       = [CEV0025; CEV002; CEV001; CEVN2005; CEVN2009];
        FRCStopCritVal(1+(i-1)*lengthStopCrits:i*lengthStopCrits)       = [FRC0025; FRC002; FRC001; FRCN2005; FRCN2009];
        M0StopCritVal(1+(i-1)*lengthStopCrits:i*lengthStopCrits)        = [M0Stop0025; M0Stop002; M0Stop001; M0StopN2005; M0StopN2009];
        M1StopCritVal(1+(i-1)*lengthStopCrits:i*lengthStopCrits)        = [M1Stop0025; M1Stop002; M1Stop001; M1StopN2005; M1StopN2009];
        M2StopCritVal(1+(i-1)*lengthStopCrits:i*lengthStopCrits)        = [M2Stop0025; M2Stop002; M2Stop001; M2StopN2005; M2StopN2009];
        LCIStopCritBreath(1+(i-1)*lengthStopCrits:i*lengthStopCrits)    = [breath0025; breath002; breath001; breathN2005; breathN2009];
        LCIStopCritRatio(1+(i-1)*lengthStopCrits:i*lengthStopCrits)     = [ratio0025; ratio002; ratio001; ratioN2005; ratioN2009];       
        
    end
    header = outputMeasurement(1,:);
    header = header(~cellfun(@isempty,header));
    
    LCIStopCrit             = cell(5, length(LCIStopCritVal));
    LCIStopCrit(1, 1:lengthStopCrits:end) = header;
    LCIStopCrit(2,:)        = repmat(stopCrits, 1, numOfMethods);
    LCIStopCrit(3,:)        = num2cell(LCIStopCritVal);   
    LCIStopCrit(4,:)        = num2cell(LCIStopCritRel);
    LCIStopCrit(5,:)        = num2cell(CEVStopCritVal);
    LCIStopCrit(6,:)        = num2cell(FRCStopCritVal);    
    LCIStopCrit(7,:)        = num2cell(M0StopCritVal);
    LCIStopCrit(8,:)        = num2cell(M1StopCritVal);
    LCIStopCrit(9,:)        = num2cell(M2StopCritVal);
    LCIStopCrit(10,:)        = num2cell(LCIStopCritBreath);   
    LCIStopCrit(11,:)        = num2cell(LCIStopCritRatio);   
end

function logic = getPossibleEndpoints(relDiff, minRelDiff)

    logic    = abs(relDiff)<minRelDiff;
    logic    = [0; 0; logic(3:end).*logic(2:end-1).*logic(1:end-2)];     %Find the first index where 3 rel differences are below minRelDiff
         
end

function logic = getValAcordingToStopCriteria(val, minRelDiff)

    relDiff1 = abs(val(1:end-2)-val(2:end-1));
    relDiff2 = abs(val(1:end-2)-val(3:end));
    
    logic    = [ 0; 0; val(1:end-2)<minRelDiff.*relDiff1<minRelDiff.*relDiff2<minRelDiff];
end

function [varargout] = nanIfEmpty(varargin)
           
    varargout = varargin;
    
    for i = 1:nargin
        if isempty(varargin{i})
            varargout{i} = nan;
        end
    end   
end