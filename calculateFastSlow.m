function gas = calculateFastSlow(table, gas, validData, parameters)
    addpath('.\FMINSEARCHBND');                          %Add the fminsearchbnd function
    graphState = parameters.Simulation.graphState;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Set variables for x, y Axis
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    xCEV = gas.General.cev_ds(validData==1)*1000;                   XAxisCEV = 'CEV [L]';
    xNum = (1:length(gas.General.cev_ds(validData==1)))';           XAxisNum = 'Breath Number';
    xTO  = gas.General.cev_ds(validData==1)./gas.General.frcao(validData==1); XAxisTO = 'TO';
    yCet = gas.General.cet_norm(validData==1);                      YAxisCet = 'N2 Concentration [%]';
    yVN2 = (gas.exp(validData==1)-gas.reInsp(validData==1))*1000;   YAxisVN2 = ['Expired ' gas.name ' Volume [L]'];

    indices   = find(yCet<0.025);
    targetIndexConcentration=getLCIIndex(indices,0);
    
    x0Cet  = [0.8; -0.615; 0.2; -0.4; 0.9];            %Initial search parameters for concentration measurements
    x0VN2  = [0.8; -0.75; 0.2; -0.25; 0.9];            %Initial search parameters for volume measurements
    
    if graphState == 1
        %h1 = getOrMakeFigure('Fast Slow Log functions Exponential N2Cet/#');
        h2 = getOrMakeFigure('Fast Slow Log functions Exponential N2CetNorm/CEV');
    else 
        h1 = 0;
        h2 = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fast slow fit for whole measurement N2[%]/#
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %axisPair='Exponential N2Cet/#';
    %expFuncitonsC   = calculateFastSlowFunction(gas, xNum, yCet, XAxisNum, YAxisCet, x0Cet,  axisPair, h1, graphState, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fast slow fit for whole measurement VolExpN2/#
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %axisPair='Exponential VolExpN2/#';
    %expFuncitonsV   = calculateFastSlowFunction(gas, xNum, yVN2, XAxisNum, YAxisVN2, x0VN2, axisPair, h2, graphState, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fast slow fit for whole measurement VolExpN2/CEV
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %axisPair='Exponential N2Cet/CEV';
    %expFuncitonsCC   = calculateFastSlowFunction(gas, xCEV, yCet, XAxisCEV, YAxisCet, repmat(x0Cet, 1, 7), axisPair, h2, graphState, 1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Fast slow fit for whole measurement VolExpN2/TO
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    axisPair='Exponential N2Cet/CEV';
    expFuncitonsCC   = calculateFastSlowFunction(gas, xTO, yCet, XAxisTO, YAxisCet, repmat(x0Cet, 1, 7), axisPair, h2, graphState, 1);
    
    
    %expFunctions = [expFuncitonsV'; expFuncitonsC'; expFuncitonsCC'];
    expFunctions = [expFuncitonsCC'];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %  Calculate further parameters from the functions found
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    expFVol = cell(length(expFunctions),1); 
    for i = 1:length(expFuncitonsCC)
        %expFVol{1+3*(i-1)} = fastSlowMBWexp(gas, xNum, yCet, expFuncitonsC(i), parameters, targetIndexConcentration);
        %expFVol{2+3*(i-1)} = fastSlowMBWexp(gas, xNum, yVN2, expFuncitonsV(i), parameters, targetIndexConcentration);
        expFVol{i} = fastSlowMBWexp(gas, xTO, yCet, expFuncitonsCC(i), parameters, targetIndexConcentration);
    end
    gas.fastSlow = expFVol;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Prediction part
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if parameters.Simulation.fastSlowPredict == 1 && gas.Standard.lci(1) ~= 0       %Only calculate the prediction if a lci 2.5% was found else no relative error and so on can be calculated
        nOfBreathsUsed   = 8:length(xNum);   %Get breaths from all to minimally 5
        predictionStructs = cell(length(nOfBreathsUsed),1);
        x0CetCEV = repmat(x0Cet, 1, 7);  
        wh = waitbar(0,'Calculating Washout Predictions', 'Name','Predicting Washout');
        for i = nOfBreathsUsed
            waitbar(i/length(nOfBreathsUsed));  
            [predictionStructs{i-nOfBreathsUsed(1)+1,1},  x0CetCEV] = predictLCI(gas, xCEV(1:i), yCet(1:i), XAxisCEV, YAxisCet, x0CetCEV, 'Concentration/CEV',           h2, graphState);
        end
        close(wh);
        
        %New version of the prediction, working yet and quite simplyfied
        [ outputModel, outputModelInterp, outputParams, LCIStopCrit] = calculateAndPreparePredictionDataV2( parameters, table, gas, predictionStructs, validData );
        savePredictionDataV2( parameters, table, gas, outputModel, outputModelInterp, outputParams, LCIStopCrit);
    end
end

%
% Input: gas:       a Gas Object containing the fields gas.General.cev_ds and
%                   gas.General.frcao 
%        validData: Breathsa which shall be included in LCI calculation
%        x:         A vector containing the x values used to describe the
%                   fast slow function 
%        y:         A vector containing the y values used to describe the
%                   fast slow function 
%        XAxisString & yAxisString: Strings to annotate the plot       
%        x0:        A vecot containing the starting point of the fminsearch 
%                   function
%        method:    A string to annotate the plots, which method is used
%        h:         A figure handle to display the fast/slow fits
%        graphStat: If graphs shall be plotted
%        version:   Fast/Slow(1) or Prediction(2)
%
function structFunctions = calculateFastSlowFunction(gas, x, y, XAxisString, YAxisString, x0, axisPair, h, graphState, version)    
       
    exp1   = @(a,x)(a(1).*exp(a(2)*x));
    exp2   = @(a,x)(a(1).*exp(a(2)*x)+a(3).*exp(a(4)*x));
    exp3   = @(a,x)((a(1)-a(3)).*exp(a(2)*x)+a(3).*exp(a(4)*x));
    expQ   = @(a,x)(a(1).*exp(a(2)*x.^a(3)));
  
    maxNumOfBreaths = max(150, length(x)*1.5); %Assumption when the test would have to be finished at the latest to make sure the fitted functions converge below 1/40th of the starting concentration at least until then
    maxCEV          = max(100, x(end)*1.5);
    
    
    LBexp  = [0; -20];
    UBexp  = [10; 0 ];
    LBquad = [0; -20; 0];
    UBquad = [10; 0;  1];
    
    methodCounter = 1;
    
    if graphState == 1
        set(0, 'currentfigure', h);
    end
    
    opts = optimset('MaxFunEvals',50000, 'MaxIter',10000, 'Display','off');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method Sweden
    %
    % Paper:    Slow and fast lung compartments in cystic fibrosis measured 
    %           by nitrogen multiple-breath washout
    % Appendix: Determination  of  the  fast  and  slowly  ventilated  lung
    %           compartments  
    % Link:     http://jap.physiology.org/content/117/7/720.abstract
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     if version == 1
%         firstBreath=ceil(length(x)*0.1);
%         lastBreath=min(ceil(length(x)*0.9), max([length(x)-4, firstBreath]));
%         rSqr=zeros(size(firstBreath:lastBreath));
% 
%         ffinds = cell(length(firstBreath:lastBreath), 1);           
%         for i = firstBreath:lastBreath %Use at least 5 breaths to fit the curve       
%             OLS = @(b) sum(((exp1(b,x(i:end))-y(i:end))./y(i:end)).^2+(exp1(b, maxNumOfBreaths)>y(1)*0.025)*10);
%             B   = fminsearchbnd(OLS,x0(3:4),LBexp,UBexp,opts);
% 
%             rSqr(i-firstBreath+1) = getR2(x(firstBreath:end),y(firstBreath:end),exp1, B); 
%             ffinds{i-firstBreath+1} = B;
%         end
%         [~, i] = max(rSqr);         %Find best fit
%         pslow=ffinds{i};            %Get the corresponding function
% 
%         OLS = @(b) sum(((exp2(b,x)-y)./y).^2+(exp2(b, maxNumOfBreaths)>y(1)*0.025)*10);
%         LBSweden = [LBexp(1:2); pslow];
%         UBSweden = [UBexp(1:2); pslow];
%         pfast   = fminsearchbnd(OLS,[x0(1:2); pslow],LBSweden,UBSweden,opts);
% 
%         pfast = sortParams(pfast);
% 
%         %r2fast  = getR2y(y,exp1(pfast(1:2), y-exp1(pslow,y)));
%         r2      = getR2(x,y,exp2, pfast); 
%         r2first = getR2(x(1:5),y(1:5),exp2, pfast);
%         r2last  = getR2(x(end-4:end),y(end-4:end),exp2, pfast);
% 
%         structFunctions(methodCounter) = constructFastSlowFunctions(['Sweden ' axisPair], ...
%                 exp2, pfast, exp1, pfast(1:2), exp1, pslow, 0, 0, r2);
%         structFunctions(methodCounter).r2 = [r2, r2first, r2last];
%         if graphState == 1
%             h = subplot(3,3,methodCounter);
%             annotateSubplot(h, 'Method Sweden', gas, XAxisString, YAxisString, x, y, exp2,  pfast, exp1, pfast(1:2), exp1, pfast(3:4), r2);
%         end
%         methodCounter = methodCounter+1;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method direct y = (a-c)*exp(b*x)+c*exp(d*x)    
    % A function with two power terms is fitted to the data with the
    % Matlab function fit()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LB = [0.6; -0.922; 0.025; -0.527];
    UB = [1.1; -0.527; 0.6; -0.05];
    x0 = [0.8; -0.615; 0.2; -0.4];
    OLS = @(b) sum(((exp3(b,x)-y)./y).^2+(exp3(b, maxCEV)>0.025)*10);
    B   = fminsearchbnd(OLS,x0,LB,UB,opts);
    B   = sortParams(B);   
    r2  = getR2(x,y,exp3, B); 
    r2first = getR2(x(1:5),y(1:5),exp3, B);
    r2last  = getR2(x(end-4:end),y(end-4:end),exp3, B);
    
    structFunctions(methodCounter) = constructFastSlowFunctions(['Direct y = a*exp(b*x)+c*exp(d*x) ' axisPair], ...
                exp3, B, exp1, B(1:2), exp1, B(3:4), 0, 0, r2);  
    structFunctions(methodCounter).r2 = [r2, r2first, r2last];
       
    if graphState == 1
        h = subplot(3,3,methodCounter);
        annotateSubplot(h, 'Matlab Fit y = a*exp(b*x)+c*exp(d*x)', gas, XAxisString, YAxisString, x, y, exp3, B, exp1, [B(1)-B(3) B(2)], exp1, [B(3) B(4)], r2);
    end 
    methodCounter = methodCounter + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method direct y = a*exp(b*x^c)
    %
    % A function with two power terms is fitted to the data with the
    % Matlab function fit()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LB = [0.8; -2; 0.2];
    UB = [2; -0.1; 1];
    x0 = [1; -0.7; 0.8];
    OLS = @(b) sum(((expQ(b,x)-y)./y).^2+(expQ(b, maxCEV)>0.025)*10);
    B   = fminsearchbnd(OLS,x0,LB,UB,opts);
    r2  = getR2(x,y,expQ, B); 
    r2first = getR2(x(1:5),y(1:5),expQ, B);
    r2last  = getR2(x(end-4:end),y(end-4:end),expQ, B);
    structFunctions(methodCounter) = constructFastSlowFunctions(['Direct y = a*exp(b*x^c)' axisPair], ...
                expQ, B, expQ, B, @(x, p) 0, [0 0], 0, 0, r2);  
    structFunctions(methodCounter).r2 = [r2, r2first, r2last];

    if graphState == 1
        h = subplot(3,3,methodCounter);
        annotateSubplot(h, 'Matlab Fit y = a*exp(b*x^c)', gas, XAxisString, YAxisString, x, y, expQ, B, expQ, B, @(x, p)0, 0, r2);
    end   
    methodCounter = methodCounter + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method direct y = a*exp(b*x)+c*exp(d*x^f)
    %
    % A function with two power terms is fitted to the data with the
    % Matlab function fit()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LB = [0.1; -0.671; 0.1; -1.0; 0.3];
    UB = [1.0; -0.527; 1.3; -0.01; 1.0];
    x0 = [0.9; -0.6; 0.1; -0.2; 0.8];
    OLS = @(b) sum((((exp1(b(1:2),x)+expQ(b(3:5),x))-y)./y).^2 ...
            +exp1(b(1:2), x(end))+expQ(b(3:5), maxCEV)>0.025)*10;
    B   = fminsearchbnd(OLS,x0,LB,UB,opts);
    r2  = getR2(x,y,@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B); 
    r2first = getR2(x(1:5),y(1:5),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B);
    r2last  = getR2(x(end-4:end),y(end-4:end),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B);

    structFunctions(methodCounter) = constructFastSlowFunctions(['Direct y = a*exp(b*x)+c*exp(d*x^f) ' axisPair], ...
                @(b,x) exp1(b(1:2),x)+expQ(b(3:5), x), B, exp1, B(1:2), expQ, B(3:5), 0, 0, r2);  
    structFunctions(methodCounter).r2 = [r2, r2first, r2last];

    if graphState == 1
        h = subplot(3,3,methodCounter);
        annotateSubplot(h, 'Matlab Fit y = a*exp(b*x)+c*exp(d*x^f)', gas, XAxisString, YAxisString, x, y, @(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B, exp1, B(1:2), expQ, B(3:5), r2);
    end
    methodCounter = methodCounter + 1;
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method split at 45° angle
    %
    % The fit from the direct method is taken and derived. Where tangent of
    % the funcition is 45° the split between the fast and slow compartment
    % is assumend. The data prior to the 45° split is taken for the fast
    % compartment estimation, the data after the split for the slow
    % compartment contribution
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
%     index = getAngleIndex(x, y, x0(1:4),[LBexp; LBexp],[UBexp; UBexp], opts, exp2);
%     index = max(index, 5);
%     index45 = {axisPair, index};
%      
%     %Split at 45° with y = a*exp(b*x)+c*exp(d*x), calculate the slow compartment first 
%     OLS	= @(b) sum(((exp1(b,x(index+1:end))-y(index+1:end))./y(index+1:end)).^2+(exp1(b, maxNumOfBreaths)>y(1)*0.025)*10);
%     pslow  = fminsearchbnd(OLS,x0(3:4),LBexp,[max(y); 0],opts);
%     
%     OLS = @(b) sum(((exp1(b,x)-(y-exp1(pslow,x)))./(y-exp1(pslow,x))).^2+(exp1(b, maxNumOfBreaths)>y(1)*0.025)*10);
%     pfast  = fminsearchbnd(OLS,x0(1:2),LBexp,UBexp,opts);
%     
%     r2slow = getR2(x(index+1:end),y(index+1:end),exp1,pslow); 
%     r2fast = getR2(x(1:index),y(1:index),exp1, pfast);        %Calculate the R-Square value of the summed functions
%     r2  = getR2y(y, exp1(pslow, x)+exp1(pfast, x));
%     r2first  = getR2y(y(1:5), exp1(pslow, x(1:5))+exp1(pfast, x(1:5)));
%     r2last  = getR2y(y(end-4:end), exp1(pfast, x(end-4:end)));
% 
%     B = sortParams([pslow; pfast]);
%     pfast = B(1:2);
%     pslow = B(3:4);
%     structFunctions(methodCounter) = constructFastSlowFunctions(['45° slow first exp' axisPair], ...
%                                                     exp2, [pfast; pslow], exp1, pfast, exp1, pslow, 0, 0, r2); 
%     structFunctions(methodCounter).r2 = [r2, r2first, r2last];
%      
%     if graphState == 1
%         h = subplot(3,3,methodCounter);
%         annotateSubplot(h, ['45° slow first exp' axisPair], gas, XAxisString, YAxisString, x, y, exp2, [pfast; pslow], exp1, pfast, exp1, pslow, r2);  
%         hold on
%         line([x(index) x(index)], [min(y) max(y)], 'LineWidth', 1, 'Color', 'black', 'LineStyle','--'); %Mark the 45° line
%         hold off
%     end
%     methodCounter = methodCounter + 1;
%         
%     
    %%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Split at 45° with y = a*exp(b*x)+c*exp(d*x^f), calculate the slow compartment first 
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     OLS	= @(b) sum(((expQ(b,x(index+1:end))-y(index+1:end))./y(index+1:end)).^2+(expQ(b, maxNumOfBreaths)>y(1)*0.025)*10);
%     pslow  = fminsearchbnd(OLS,x0(3:5),LBquad,[max(y); 0; 1],opts);
% 
%     OLS = @(b) sum(((exp1(b,x)-(y-expQ(pslow,x)))./(y-expQ(pslow,x))).^2);
%     pfast  = fminsearchbnd(OLS,x0(1:2),LBexp,UBexp,opts);
% 
%     r2slow = getR2(x(index+1:end),y(index+1:end),expQ,pslow); 
%     r2fast = getR2(x(1:index),y(1:index),exp1, pfast);        %Calculate the R-Square value of the summed functions
%     r2  = getR2y(y, [exp1(pfast, x) + expQ(pslow, x)]);
%     r2first  = getR2y(y(1:5), exp1(pfast, x(1:5)));
%     r2last  = getR2y(y(end-4:end), expQ(pslow, x(end-4:end)));
% 
%     cutoff = x(index);
% 
%     structFunctions(methodCounter) = constructFastSlowFunctions(['45° slow first quad' axisPair], ...
%                                                     @(b,x) exp1(b(1:2),x)+expQ(b(3:5), x), [pfast; pslow], exp1, pfast, expQ, pslow, 0, 0, r2); 
%     structFunctions(methodCounter).r2 = [r2, r2first, r2last];
% 
%     if graphState == 1
%         h = subplot(3,3,methodCounter);
%         annotateSubplot(h, ['45° slow first quad' axisPair], gas, XAxisString, YAxisString, x, y,  @(b,x) exp1(b(1:2),x)+expQ(b(3:5), x), [pfast; pslow], exp1, pfast, expQ, pslow, r2);  
%         hold on
%         line([x(index) x(index)], [min(y) max(y)], 'LineWidth', 1, 'Color', 'black', 'LineStyle','--'); %Mark the 45° line
%         hold off
%     end
%     methodCounter = methodCounter + 1;
%         
     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method start fit at extreme points with constrained adaptation in the
    % fit function
    %
    % Under the assumption that the first and last datapoints of the
    % washout fit the fast and slow washout the strongest, a fit routine is
    % started with exponential functions fitted to the first datapoints.
    % The single exponential function is taken and its parameters are used
    % to find a double exponential functions with parameters close in the
    % range of the first function. 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if version == 1
%         boundaryLength=max(round(length(x)*0.33), 5);
%         bnd = 1.5;   %The parameters can be searched within a limit of bnd
% 
%         OLS     = @(b) sum(((exp1(b,x(1:boundaryLength))-y(1:boundaryLength))./y(1:boundaryLength)).^2);
%         pfast   = fminsearchbnd(OLS,x0(1:2),LBexp,UBexp,opts);
%         pfast(1) = pfast(1)*0.8;
% 
%         OLS     = @(b) sum(((exp1(b,x(end-boundaryLength:end))-y(end-boundaryLength:end))./y(end-boundaryLength:end)).^2);
%         pslow   = fminsearchbnd(OLS,x0(3:4),LBexp,UBexp,opts);
% 
%         OLS     = @(b) sum((((exp1(b(1:2),x)+expQ(b(3:5),x))-y)./y).^2+((exp1(b(1:2),maxNumOfBreaths)+expQ(b(3:5),maxNumOfBreaths))>y(1)*0.025)*10);
% 
%         LB      = [pfast(1)/bnd; pfast(2)*bnd; pslow(1)/bnd; pslow(2)*bnd; 0];
%         UB      = [pfast(1)*bnd; pfast(2)/bnd; pslow(1)*bnd; pslow(2)/bnd; 1];
%         pMixed  = fminsearchbnd(OLS,[pfast; pslow; 0.9],LB,UB,opts);
% 
%         pMixed(1:4) = sortParams(pMixed(1:4));
% 
%         r2fast  = getR2(x(1:boundaryLength),y(1:boundaryLength),exp1, pMixed(1:2)); 
%         r2slow  = getR2(x(end-boundaryLength:end),y(end-boundaryLength:end),expQ, pMixed(3:5));
%         r2      = getR2(x,y,@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), pMixed);
%         r2first  = getR2(x(1:5),y(1:5),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), pMixed);
%         r2last   = getR2(x(end-4:end),y(end-4:end),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), pMixed);
% 
%         structFunctions(methodCounter) = constructFastSlowFunctions(['Constrained ' axisPair], ...
%             @(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), pMixed, exp1, pMixed(1:2), expQ, pMixed(3:5), 0, 0, r2);
%         structFunctions(methodCounter).r2 = [r2, r2first, r2last];
% 
%         if graphState == 1
%             h = subplot(3,3,methodCounter);
%             %semilogy(x,y, '.', x, f1(pfast, x), x, f1(pslow, x), x, f2(pMixed, x));
% 
%             annotateSubplot(h, 'Constrained boundary data for fit params', gas, XAxisString, YAxisString, x, y, @(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), pMixed, exp1, pMixed(1:2), expQ, pMixed(3:5), r2);    
%         end
%         methodCounter = methodCounter + 1;
%     end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method offset
    %
    % The data has a slight offset from the x and y axis, to counteract
    % the offset the funciton is chaged from 
    % y= a*x^x+c*x^d            to 
    % y= a*x^b+c*x^d + f + g*x
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if version == 1
%         startIndex = ceil(length(x)/2);
%         OLS     = @(b) sum(((exp1(b,x(startIndex:end))-y(startIndex:end))./y(startIndex:end)).^2);
%         pslow   = fminsearchbnd(OLS,x0(3:4),LBexp,UBexp,opts);
% 
%         fOffset = @(a,x)(exp2(a(1:4),x)+a(5)*x+a(6));
% 
%         OLS = @(b) sum(((fOffset(b,x)-y)./y).^2+(fOffset(b, maxNumOfBreaths)>y(1)*0.025)*10);
%         B   = fminsearchbnd(OLS,[x0(1:2); pslow; 4e-05; 1.66e-03],[LBexp; LBexp; 0; 1e-6],[UBexp; UBexp; 1e-3; 10e-3],opts);
% 
%         B(1:4) = sortParams(B(1:4));
% 
%         r2      = getR2(x,y,fOffset, B);
%         r2first = getR2(x(1:5),y(1:5),fOffset, B);
%         r2last  = getR2(x(end-4:end),y(end-4:end),fOffset, B);
% 
%         structFunctions(methodCounter) = constructFastSlowFunctions(['Offset ' axisPair], ...
%                     fOffset, B, exp1, B(1:2), exp1, B(3:4), @(a,x) a(1)*x+a(2), B(5:6), r2);
%         structFunctions(methodCounter).r2 = [r2, r2first, r2last];
%         
%         if graphState == 1
%             h = subplot(3,3,methodCounter);
%             %semilogy(x,y, '.', x, f1(B(1:2),x), x, f1(B(3:4),x), x, fOffset(B,x));
% 
%             annotateSubplot(h, ['Method Offset' axisPair], gas, XAxisString, YAxisString, x, y, exp2, B, exp1, B(1:2), exp1, B(3:4), r2);
%             text(min(x)*1.2,max(y)*0.55, ['a: ' num2str(B(5)) ''])
%             text(min(x)*1.2,max(y)*0.35, ['b: ' num2str(B(6)) ''])        
%         end
%         methodCounter = methodCounter + 1;
%     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method FRC dependent equal amount
    %
    % Under the assumption that the contribution of the fast and slow
    % compartments are the same at LCI/2 the best fitting function is
    % found with the boundary conditions:
    % a*e^(b*(xLCI)) -c*e^(d*(xLCI)) = 0
    % a*e^(b*x) - c*e^(d*x))         = y
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     if version == 1
%         LCI = gas.General.cev_ds./gas.General.frcao;
%         %LCI(validData==0) = [];
%         index = find(LCI>max(LCI)/2, 1);                        %find the index where the LCI is over half the maximum, and the fast and slow funciton should intersect
%         index = max(index, 5); 
%         index = min(index, max(length(x)-5,2));
% 
%         OLS     = @(b) sum(((exp1(b,x(index+1:end))-y(index+1:end))./y(index+1:end)).^2);
%         pslow   = fminsearchbnd(OLS,B(3:4),LBexp,UBexp,opts);
% 
%         LCI = mean(x(index-1:index+1));
%         weightOfContributionEquality=length(x); 
% 
%         OLS = @(b) sum((((exp2(b,x)+ ...
%                         weightOfContributionEquality*(exp1(b(1:2), LCI) - exp1(b(3:4), LCI)))-y)./y).^2+(exp2(b, maxNumOfBreaths)>y(1)*0.025)*10);           
%         pEqual   = fminsearchbnd(OLS,[x0(1:2); pslow],[LBexp; LBexp],[UBexp; UBexp],opts);
% 
%         pEqual = sortParams(pEqual);
% 
%         r2       = getR2(x,y,exp2, pEqual); 
%         r2first  = getR2(x(1:5),y(1:5),exp2, pEqual);
%         r2last   = getR2(x(end-4:end),y(end-4:end),exp2, pEqual);
% 
%         structFunctions(methodCounter) = constructFastSlowFunctions(['Equal ' axisPair], ...
%                     exp2, pEqual, exp1, pEqual(1:2), exp1, pEqual(3:4), 0, 0, r2);
%         structFunctions(methodCounter).r2 = [r2, r2first, r2last];
%         
%         if graphState == 1
%             h = subplot(3,3,methodCounter);
%             %semilogy(x,y, '.', x, f1(pEqual(1:2), x), x, f1(pEqual(3:4), x), x, f2(pEqual, x));
% 
%             annotateSubplot(h, 'Equal contribution at LCI<min(LCI)*2', gas, XAxisString, YAxisString, x, y, exp2, pEqual, exp1, pEqual(1:2), exp1, pEqual(3:4), r2);          
%             hold on
%             line([x(index) x(index)], [min(y) max(y)], 'LineWidth', 1, 'Color', 'black', 'LineStyle','--');
%             hold off
%         end
%         methodCounter = methodCounter + 1;
%     end    
end

%Get R2 based on function values
function r2 = getR2(x,y,f, p)
    r2     = 1-sum((y-f(p,x)).^2)/sum((y-mean(y)).^2);
end

%Get R2 based on measured y and calculated y (yf)
function r2 = getR2y(y, yf)
    r2     = 1-sum((y-yf).^2)/sum((y-mean(y)).^2);
end

%Get R2 with relative errors
function r2 = getR2CV(x,y,f, p)
    r2     = 1-sum(((y-f(p,x)).^2)./y)/sum(((y-mean(y)).^2)./y);
end

%Get R2 with relative errors based on measured y and calculated y (yf)
function r2 = getR2yCV(y, yf)
    r2     = 1-sum(((y-yf).^2)./y)/sum(((y-mean(y)).^2)./y);
end

function index = getAngleIndex(x, y, x0, LB, UB, opts, exp2)
    OLS = @(b) sum(((exp2(b,x)-y)./y).^2);
    B   = fminsearchbnd(OLS,x0,LB,UB,opts);
    
    yfit=exp2(B,x);
    derivNorm=diff(yfit./max(yfit)*length(yfit));%./diff(x);          %Derivate the normalized fit function
    angle=45;
    angle=cot(angle/180*pi());                                      %Set the slope at which the fitted function shall be splitted 
    %[~, index] = min(abs(derivNorm+angle));                         %Find the 45° angle of a strictly monotonuous falling function, and set the split location
    index   = find(diff(sign(derivNorm+angle)), 1)+1;
    if isempty(index)
        index = 0;
        warning('No 45° Angle found')
    end
    index   = min(max(index, 5), length(x)-5);                   %Set the meeting point at least 5 breaths from the start or end
end

%With f1, p1, x01, LB1, UB1 for the fast compartment and all parameters
%with index 2 for the slow compartment
function [p1, p2, r21, r22, r2] = getAngleParametersSlowCompartmentFirst(x, y, index, f1, f2, x01, x02, LB1, LB2, UB1, UB2, opts)
    
    OLS	= @(b) sum(((f1(b,x(index+1:end))-y(index+1:end))./y(index+1:end)).^2);
    p1  = fminsearchbnd(OLS,x01,LB1,UB1,opts);
    
    OLS = @(b) sum(((f2(b,x(1:index))-(y(1:index)-f1(p1,x(1:index))))./(y(1:index)-f1(p1,x(1:index)))).^2);
    p2  = fminsearchbnd(OLS,x02,LB2,UB2,opts);
    
    r21 = getR2(x(index+1:end),y(index+1:end),f1,p1); 
    r22 = getR2(x(1:index),y(1:index),f2, p2);        %Calculate the R-Square value of the summed functions
    r2  = getR2y(y, [f1(p1, x(1:index))+f2(p2, x(1:index));f2(p2, x(index+1:end))]);
end

%With f1, p1, x01, LB1, UB1 for the slow compartment and all parameters
%with index 2 for the fast compartment
function [p1, p2, r21, r22, r2] = getAngleParametersFastCompartmentFirst(x, y, index, f1, f2, x01, x02, LB1, LB2, UB1, UB2, opts)
    
    OLS	= @(b) sum(((f1(b,x(1:index))-y(1:index))./y(1:index)).^2);
    p1  = fminsearchbnd(OLS,x01,LB1,UB1,opts);
    
    OLS = @(b) sum(((f2(b,x(index+1:end))-(y(index+1:end)-f1(p1,x(index+1:end))))./(y(index+1:end)-f1(p1,x(index+1:end)))).^2);
    p2  = fminsearchbnd(OLS,x02,LB2,UB2,opts);
    
    r21 = getR2(x(1:index),y(1:index),f1,p1); 
    r22 = getR2(x(index+1:end),y(index+1:end),f2, p2);        %Calculate the R-Square value of the summed functions
    r2  = getR2y(y, [f1(p1, x(1:index))+f2(p2, x(1:index));f2(p2, x(index+1:end))]);
end

%
%functtion to display a subplot of a fast slow decomposed LCI
%
%input: h           = handle to a subplot
%       titleS      = title of the subplot
%       x           = x axis data
%       y           = y axis data
%       a, b, c, d  = fitted function a*e^b+c*e^d whereby the first
%                     exponential is for the slow compartment and the
%                     second one for the fast
%       r2          = r squared, showing how good the fit matches the data
%
function annotateSubplot(h, titleS, gas, XAxisString, YAxisString, x, y, f, p, ff, fp, sf, sp, r2)
    subplot(h)
    semilogy(x,y, 'b.', x, sf(sp,x), 'g-', x, ff(fp,x), 'r-', x, f(p,x), 'b-');
    title(titleS);
    axis([0, max(x)*1.05, min(y)*0.5, max(y)*1.05]) 
    legend('Breaths', 'Slow Compartment', 'Fast Compartment', 'Fitted Curve');
    xlabel(XAxisString);
    ylabel(YAxisString); 
    text(min(x)*1.2,max(y)*0.7, ['R-Squared: ' num2str(r2) '']);
    if gas.Standard.nBreaths(1) <= length(x)
        hold on
        semilogy(x(gas.Standard.nBreaths(1)), y(gas.Standard.nBreaths(1)), 'ko');
        hold off
    end
end

%Output:    Function
%           FRC, V fast and slow 
           
function data = fastSlowMBWexp(gas, x, y, f, parameters, targetIndices)
    
    DSpre     =   parameters.Device.volumePrecap;
    

    paramFast = f.pFast;
    paramSlow = f.pSlow;
    
    ffast = @(x)f.fFast(paramFast, x);
    fslow = @(x)f.fSlow(paramSlow, x);
            
    %Get LCI 2.5% Index
    
    if targetIndices == 0   %No 2.5% was reached
        VGasSlow= 0;
        VGasFast= 0;
        FRCFast  = 0;
        FRCSlow  = 0;
        w     = 0;
        wFast = 0;
        wSlow = 0;
        AlveolarTidalVolumeFast = 0;
        AlveolarTidalVolumeSlow = 0;
        RegionalSpecificVentilationFast = 0;
        RegionalSpecificVentilationSlow = 0;
        FunctionalDeadSpaceVolume     = 0;
        FunctionalDeadSpaceVolumeFast = 0;        
        FunctionalDeadSpaceVolumeSlow = 0;
        meanWhole   = 0;
        meanFirst   = 0;
        meanLast    = 0;
        meanMiddle  = 0;
        stdWhole    = 0;
        stdFirst    = 0;
        stdLast     = 0;
        stdMiddle   = 0;       
        resFrontToBack = 0;            
        signWhole   = 0;
        signFirst   = 0;
        signLast    = 0;
        signMiddle  = 0;
        residuals   = 0;
        
        warning('No LCI 2.5% found so no valid calculation for the fast slow decomposition was made');
    else
        
        fastSlowRatioAtTarget = ffast(x(targetIndices)) / (fslow(x(targetIndices)) + ffast(x(targetIndices)));
    
        %Fitted Gas volume function
        fVGasFast = @(x)gas.General.netExpVol(x)'.*ffast(x)./(ffast(x)+fslow(x));
        fVGasSlow = @(x)gas.General.netExpVol(x)'.*fslow(x)./(ffast(x)+fslow(x));

        %Gas Volume
        VGasTotal     = gas.General.cumNetExpVol(targetIndices);
        VGasSlow      = sum(fVGasSlow(1:targetIndices));
        VGasFast      = VGasTotal - VGasSlow;

        %Gas concentration in the beginning of the measurement and the end(2.5%)
        CetGasInitial = gas.General.cet_start(1);
        CetGasFinal   = gas.et(targetIndices);

        %FRC calculated for total, fast and slow lung compartments
        FRCTotal = gas.Standard.frc(1);
        FRCFast  = VGasFast/(CetGasInitial - CetGasFinal * fastSlowRatioAtTarget ) - DSpre * fastSlowRatioAtTarget;
        FRCSlow  = FRCTotal - FRCFast;

        %Washout coefficient
        w     = nthroot(max(gas.General.netExpVol(targetIndices)/gas.General.netExpVol(1), 0), targetIndices);
        wFast = nthroot(max(fVGasFast(targetIndices)/fVGasFast(1), 0), targetIndices);
        wSlow = nthroot(max(fVGasSlow(targetIndices)/fVGasSlow(1), 0), targetIndices); 
        %If the contribution remains the same over the whole measurement,
        %it is possible that with the calculation above w is > 1 which
        %would imply that the slow regions contribute more VN2 the longer
        %the test runs which should be impossible, therefore the formula
        %below is used in that case. 
        if wSlow > 1 
            wSlow = nthroot(fslow(targetIndices)/fslow(1), targetIndices);
        end
        %Alveolar Tidal Volume
        AlveolarTidalVolumeFast = (FRCFast/wFast)-FRCFast;      %rVtA,Speed
        AlveolarTidalVolumeSlow = (FRCSlow/wSlow)-FRCSlow;
     
        %Regional Specific Ventilation
        RegionalSpecificVentilationFast = AlveolarTidalVolumeFast/FRCFast; %rV'/VSpeed
        RegionalSpecificVentilationSlow = AlveolarTidalVolumeSlow/FRCSlow;
        
        %Total expired N2
        CeN2     = sum(gas.General.netExpVol(1:targetIndices))/gas.General.cev_ds(targetIndices);
        CeN2Fast = sum(fVGasFast(1:targetIndices))/gas.General.cev_ds(targetIndices);
        CeN2Slow = sum(fVGasSlow(1:targetIndices))/gas.General.cev_ds(targetIndices);
        
        %FunctionalDeadSpace Volume 
        FunctionalDeadSpaceVolume       = (1 - CeN2/(polyval([ones(1, targetIndices) 0], w)/targetIndices));
        FunctionalDeadSpaceVolumeFast   = (1 - CeN2Fast/(polyval([ones(1, targetIndices) 0], wFast)/targetIndices));
        FunctionalDeadSpaceVolumeSlow   = (1 - CeN2Slow/(polyval([ones(1, targetIndices) 0], wSlow)/targetIndices));
        
        %Calculate Residuals
        fitFunction = @(x)ffast(x)+fslow(x); 
        residuals   = (y-fitFunction(x))./y;                     %relative error calculation
        
        lr          = length(residuals);
        
        first5      = 1:5;
        last5       = lr-9:lr;
        middle5     = round(lr/2)-2:round(lr/2)+2;
        
        meanWhole   = mean(residuals);
        meanFirst   = mean(residuals(first5));
        meanLast    = mean(residuals(last5));
        meanMiddle  = mean(residuals(middle5));
        
        stdWhole    = std(residuals);
        stdFirst    = std(residuals(first5));
        stdLast     = std(residuals(last5));
        stdMiddle   = std(residuals(middle5));
        
        resFrontToBack = mean(abs(residuals(1:floor(lr/2))))/mean(abs(residuals(ceil(lr/2):end)));     %Compare the residuals from the first half vs the last half
        
        signWhole   = sum(sign(residuals))/lr;
        signFirst   = sum(sign(residuals(first5)));
        signLast    = sum(sign(residuals(last5)));
        signMiddle  = sum(sign(residuals(middle5)));
                
    end
    
    
    %Set Output structure
    data.MethodName =  f.name;
    data.FitFunction=  f.function;
            
    for i = 1:length(f.params)
        data.(['param_' num2str(i)'']) = f.params(i);
    end
                                                
    data.r2             = f.r2(1); 
    data.r2first        = f.r2(2); 
    data.r2last         = f.r2(3); 
    data.GasVolumeSlow  = VGasSlow*1e3; 
    data.GasVolumeFast  = VGasFast*1e3;   
    data.FRCSlow        = FRCSlow*1e3; 
    data.FRCFast        = FRCFast*1e3;    
    data.W              = w;                                              
    data.WSlow          = wSlow;        
    data.WFast          = wFast;          
    data.AlveolarTidalVolumeSlow        = AlveolarTidalVolumeSlow*1e3;      
    data.AlveolarTidalVolumeFast        = AlveolarTidalVolumeFast*1e3;      
    data.ReginalSpecificVentilationSlow = RegionalSpecificVentilationSlow; 
    data.ReginalSpecificVentilationFast = RegionalSpecificVentilationFast;  
    data.FunctionalDeadSpaceVolume      = FunctionalDeadSpaceVolume;        
    data.FunctionalDeadSpaceVolumeFast  = FunctionalDeadSpaceVolumeFast;    
    data.FunctionalDeadSpaceVolumeSlow  = FunctionalDeadSpaceVolumeSlow;    
    data.MeanFrontToBack                = resFrontToBack;                   
    data.MeanResiduals        = meanWhole;    data.StdResiduals        = stdWhole;  data.SignResiduals        = signWhole;  
    data.MeanFirst5Residuals  = meanFirst;    data.StdFirst5Residuals  = stdFirst;  data.SignFirst5Residuals  =  signFirst;  
    data.MeanLast10Residuals  = meanLast;     data.StdLast10Residuals  = stdLast;   data.SignLast10Residuals  =  signLast;   
    data.MeanMiddle5Residuals = meanMiddle;   data.StdMiddle5Residuals = stdMiddle; data.SignMiddle5Residuals =  signMiddle;
    
    maxResiduals = 50;
    
    residuals = padarray(residuals(1:min(length(residuals), maxResiduals)), [0, max(maxResiduals-length(residuals), 0)], 'post');
    residuals(residuals==0) = nan;
    
    for i = 1:maxResiduals
       data.(['Residual_' num2str(i, '%02d')]) = residuals(i);
    end
end
%returns fast first
function pf2 = sortParams(pf2)
    if pf2(2)>pf2(4)
       pf2=[pf2(3:4); pf2(1:2)]; 
    end
end

function s = constructFastSlowFunctions(name, fHandle, fParams, fastHandle, fastParams, slowHandle, slowParams, remHandle, remParams, r2)
   
   s = struct('name',      name,        ...
              'function',  fHandle,     ...
              'params',    fParams,     ...
              'fFast',     fastHandle,  ...
              'pFast',     fastParams,  ...
              'fSlow',     slowHandle,  ...
              'pSlow',     slowParams,  ...
              'fRem',      remHandle,   ...
              'pRem',      remParams,   ...
              'r2',        r2           );
end

