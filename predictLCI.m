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
%        0.025 = exp(-0.527.*7)     %decay constant for fast compartment
%        0.025 = exp(-0.671.*5.5)   %decay constant for fast compartment
%        0.025 = exp(-0.492.*7.5)
%        0.025 = exp(-0.738.*5)
%        0.025 = exp(-0.922.*4)
%        0.025 = exp(-0.615.*6)

function [structFunctions, x0]= predictLCI(gas, x, y, XAxisString, YAxisString, x0, axisPair, h, graphState)    
       
    exp1   = @(a,x)(a(1).*exp(a(2)*x));
    exp2   = @(a,x)(a(1).*exp(a(2)*x)+a(3).*exp(a(4)*x));
    exp3   = @(a,x)((a(1)-a(3)).*exp(a(2)*x)+a(3).*exp(a(4)*x));
    expQ   = @(a,x)(a(1).*exp(a(2)*x.^a(3)));
      
    LBexp  = [0; -10];
    UBexp  = [5; 0 ];
    LBquad = [0; -10; 0];
    UBquad = [5; 0;  1.5];
    
    maxCEV = max(100, x(end)*1.5); %Assumption when the test would have to be finished at the latest to make sure the fitted functions converge below 1/40th of the starting concentration at least until then
        
    methodCounter = 1;
    
    if size(x0,2) == 1          %If the matrix given is set generally for all functions
        x0 = repmat(x0, 1,5);
    end
    
    if graphState == 1
        set(0, 'currentfigure', h);
    end
    
    opts = optimset('MaxFunEvals',50000, 'MaxIter',10000, 'Display','off');
   
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method direct y = a*exp(b*x)+c*exp(d*x)
    %
    % A function with two power terms is fitted to the data with the
    % Matlab function fit()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LB = [0.6; -1;   0.025; -0.5];
    UB = [1.1; -0.5; 0.6;   -0.05];
    OLS = @(b) sum(((exp3(b,x)-y)./y).^2+(exp3(b, maxCEV)>0.025)*10);
    B   = fminsearchbnd(OLS,x0(1:4, 1),LB,UB,opts);
    B   = sortParams(B);   
    r2  = getR2(x,y,exp3, B); 
    r2first = getR2(x(1:5),y(1:5),exp3, B);
    r2last  = getR2(x(end-4:end),y(end-4:end),exp3, B);
    
    structFunctions(methodCounter) = constructFastSlowFunctions(['Direct y = a*exp(b*x)+c*exp(d*x) ' axisPair], ...
                exp3, B, exp1, [B(1)-B(3) B(2)], exp1, B(3:4), 0, 0, r2);  
    structFunctions(methodCounter).r2 = [r2, r2first, r2last];
       
    if graphState == 1
        h = subplot(3,3,methodCounter);
        annotateSubplot(h, 'Matlab Fit y = a*exp(b*x)+c*exp(d*x)', gas, XAxisString, YAxisString, x, y, exp3, B, exp1, [B(1)-B(3) B(2)], exp1, [B(3) B(4)], r2);
    end 
    if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
        x0(:,1)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
    else
        x0(1:4, 1) = B;
    end
    methodCounter = methodCounter + 1;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method direct y = a*exp(b*x^c)
    %
    % A function with two power terms is fitted to the data with the
    % Matlab function fit()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LB = [0.9; -3; 0.2];
    UB = [2.5; -0.2; 1];
    OLS = @(b) sum(((expQ(b,x)-y)./y).^2+(expQ(b, maxCEV)>0.025)*10+(expQ(b, 0)>0.95 && expQ(b, 0)<1.05)*10);
    B   = fminsearchbnd(OLS,x0(3:5, 2),LB,UB,opts);
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
    if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
        x0(:,2)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
    else
        x0(3:5, 2) = B;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Method direct y = a*exp(b*x)+c*exp(d*x^f)
    %
    % A function with two power terms is fitted to the data with the
    % Matlab function fit()
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    LB = [0.6; -1.5; 0.025; -1.5; 0.2];
    UB = [1.1; -0.3; 2.0; -0.2; 1.0];
    
    OLS = @(b) sum((((exp1(b(1:2),x)+expQ(b(3:5),x))-y)./y).^2 ...
            +(exp1(b(1:2), x(end))+expQ(b(3:5), maxCEV)>0.025)*10 ...
            +((exp1(b(1:2), 0)+expQ(b(3:5), 0) > 0.95) && (exp1(b(1:2), 0)+expQ(b(3:5), 0)<1.05))*10);
    B   = fminsearchbnd(OLS,x0(:,3),LB,UB,opts);
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
    if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
        x0(:,3)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
    else
        x0(:,3) = B;
    end
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Method direct y = (a-b)*exp(b*x)+b*exp(d*x)
%     %
%     % A function with two power terms is fitted to the data with the
%     % Matlab function fit()
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LB = [0.95; -2.0; 0.1; -1.0];
%     UB = [1.05; -0.2; 0.5; -0.1];
%     exp2lim   = @(a,x)((a(1)-a(3)).*exp(a(2)*x)+a(3).*exp(a(4)*x));
%     OLS = @(b) sum(((exp2lim(b,x)-y)./y).^2+(exp2lim(b, maxCEV)>0.025)*10);
%     B   = fminsearchbnd(OLS,x0(1:4, 1),LB,UB,opts);
%     B   = sortParams(B);   
%     r2  = getR2(x,y,exp2, B); 
%     r2first = getR2(x(1:5),y(1:5),exp2, B);
%     r2last  = getR2(x(end-4:end),y(end-4:end),exp2, B);
%     
%     structFunctions(methodCounter) = constructFastSlowFunctions(['Direct y = (a-c)*exp(b*x)+c*exp(d*x) ' axisPair], ...
%                 exp2, B, exp1, B(1:2), exp1, B(3:4), 0, 0, r2);  
%     structFunctions(methodCounter).r2 = [r2, r2first, r2last];
%        
%     if graphState == 1
%         h = subplot(3,3,methodCounter);
%         annotateSubplot(h, 'Matlab Fit y = (a-c)*exp(b*x)+c*exp(d*x)', gas, XAxisString, YAxisString, x, y, exp2, B, exp1, B(1:2), exp1, B(3:4), r2);
%     end 
%     methodCounter = methodCounter + 1;
%     if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
%         x0(:,4)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%     else
%         x0(1:4, 4) = B;
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Method direct y = (a-c)*exp(b*x)+c*exp(d*x^f)
%     %
%     % A function with two power terms is fitted to the data with the
%     % Matlab function fit()
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LB = [0.95; -3.0; 0.1; -1.5; 0.3];
%     UB = [1.05; -0.2; 0.7; -0.2; 1.0];
%     OLS = @(b) sum((((exp1([b(1)-b(3) b(2)],x)+expQ(b(3:5),x))-y)./y).^2 ...
%             +(exp1([b(1)-b(3) b(2)], x(end))+expQ(b(3:5), maxCEV)>0.025)*10);
%     B   = fminsearchbnd(OLS,x0(:,3),LB,UB,opts);
%     B(1) = B(1)-B(3);
%     r2  = getR2(x,y,@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B); 
%     r2first = getR2(x(1:5),y(1:5),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B);
%     r2last  = getR2(x(end-4:end),y(end-4:end),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B);
% 
%     structFunctions(methodCounter) = constructFastSlowFunctions(['Direct y = (a-c)*exp(b*x)+c*exp(d*x^f) ' axisPair], ...
%                 @(b,x) exp1(b(1:2),x)+expQ(b(3:5), x), B, exp1, B(1:2), expQ, B(3:5), 0, 0, r2);  
%     structFunctions(methodCounter).r2 = [r2, r2first, r2last];
% 
%     if graphState == 1
%         h = subplot(3,3,methodCounter);
%         annotateSubplot(h, 'Matlab Fit y = (a-c)*exp(b*x)+c*exp(d*x^f)', gas, XAxisString, YAxisString, x, y, @(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B, exp1, B(1:2), expQ, B(3:5), r2);
%     end
%     methodCounter = methodCounter + 1;
%     if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
%         x0(:,5)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%     else
%         x0(:,5) = B;
%     end
%     
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % Method split at 45° angle
%     %
%     % The fit from the direct method is taken and derived. Where tangent of
%     % the funcition is 45° the split between the fast and slow compartment
%     % is assumend. The data prior to the 45° split is taken for the fast
%     % compartment estimation, the data after the split for the slow
%     % compartment contribution
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LB = [0; -2.5;   0.1; -1.0];
%     UB = [1; -0.2;   0.6; -0.1];
%     index = getAngleIndex(x, y, x0(1:4, 4),LB,UB, opts, exp2);
%     index = max(index, 5);
%      
%     %Split at 45° with y = a*exp(b*x)+c*exp(d*x), calculate the slow compartment first 
%     LB = [0.6; -2.0;   0.1; -1.0];
%     UB = [1;   -0.4;   0.6; -0.1];
%     OLS	= @(b) sum(((exp1(b,x(index+1:end))-y(index+1:end))./y(index+1:end)).^2+(exp1(b, minCEV)>0.025)*10);
%     pslow  = fminsearchbnd(OLS,x0(3:4, 4),LB(3:4),UB(3:4),opts);
%     
%     OLS = @(b) sum(((exp1(b,x)-(y-exp1(pslow,x)))./(y-exp1(pslow,x))).^2 ...
%                    +(exp2([b; pslow], minCEV)>0.025)*10 ...
%                    +((exp2([b; pslow], 0)>0.95) && (exp2([b; pslow], 0)>1.05))*10);
%     pfast  = fminsearchbnd(OLS,x0(1:2, 4),LB(1:2),UB(1:2),opts);
%     
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
%     if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
%         x0(:,4)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%     else
%     	x0(1:4,4) = B;    
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%
%     %
%     % Split at 45° with y = a*exp(b*x)+c*exp(d*x^f), calculate the slow compartment first 
%     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LB = [0.5; -2.5; 0.1; -1.0; 0.6];
%     UB = [1.1; -1.0; 0.5; -0.3; 1];
%     OLS	= @(b) sum(((expQ(b,x(index+1:end))-y(index+1:end))./y(index+1:end)).^2+(expQ(b, minCEV)>0.025)*10);
%     pslow  = fminsearchbnd(OLS,x0(3:5, 5),LB(3:5),UB(3:5),opts);
% 
%     OLS = @(b) sum(((exp1(b,x)-(y-expQ(pslow,x)))./(y-expQ(pslow,x))).^2 ...
%                     +(exp1(b,minCEV)+expQ(pslow, minCEV)>0.025)*10  ...
%                     +((exp1(b,0)+expQ(pslow, 0)>0.95) && (exp1(b,0)+expQ(pslow, 0)>1.05))*10);
%     pfast  = fminsearchbnd(OLS,x0(1:2, 5),LB(1:2),UB(1:2),opts);
% 
%     r2  = getR2y(y, exp1(pfast, x) + expQ(pslow, x));
%     r2first  = getR2y(y(1:5), exp1(pfast, x(1:5)));
%     r2last  = getR2y(y(end-4:end), expQ(pslow, x(end-4:end)));
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
%     if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
%         x0(:,5)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%     else
%         x0(:,5) = [pfast; pslow];
%     end
%     
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %
%     %
%     %  CALCULATE THE PREDICTION WITH VARPRO
%     %
%     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %
%     % Varpro y = a*exp(b*x)+c*exp(d*x)
%     %
%     LB = [-2; -1];
%     UB = [0;   0];
%     [alpha,c,wresid,resid_norm,y_est,Regression] = ...
%      varpro(y,ones(size(y))./y,real(x0([2 4],6)),2,@(alpha)adaex_exp(alpha,x),LB,UB,opts);
%     B   = sortParams(real([c(1); alpha(1); c(2); alpha(2)]));   
%     r2  = getR2(x,y,exp2, B); 
%     r2first = getR2(x(1:5),y(1:5),exp2, B);
%     r2last  = getR2(x(end-4:end),y(end-4:end),exp2, B);
%     
%     structFunctions(methodCounter) = constructFastSlowFunctions(['Varpro y = a*exp(b*x)+c*exp(d*x) ' axisPair], ...
%                 exp2, B, exp1, B(1:2), exp1, B(3:4), 0, 0, r2);  
%     structFunctions(methodCounter).r2 = [r2, r2first, r2last];
%        
%     if graphState == 1
%         h = subplot(3,3,methodCounter);
%         annotateSubplot(h, 'Varpro Fit y = a*exp(b*x)+c*exp(d*x)', gas, XAxisString, YAxisString, x, y, exp2, B, exp1, B(1:2), exp1, B(3:4), r2);
%     end 
%     methodCounter = methodCounter + 1;
%     if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
%         x0(:,6)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%     else
%         try
%             x0(1:4, 6) = real(B);
%         catch
%             x0(:,6)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%         end
%     end
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     %
%     % Varpro y = a*exp(b*x)+c*exp(d*x^f)
%     %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     LB = [   -2;    -1; 0.01];
%     UB = [-0.01; -0.01; 0.99];
%     try
%         [alpha,c,wresid,resid_norm,y_est,Regression] = ...
%          varpro(y,ones(size(y))./y,real(x0([2 3 5],7)),2,@(alpha)adaex_pot(alpha,x),LB,UB,opts);
% 
%         B = real([c(1); alpha(1); c(2); alpha(2); alpha(3)]);
% 
%         r2  = getR2(x,y,@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B); 
%         r2first = getR2(x(1:5),y(1:5),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B);
%         r2last  = getR2(x(end-4:end),y(end-4:end),@(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B);
%     catch
%        B = [zeros(5,1)];
%        r2 = 0;
%        r2first = 0;
%        r2last = 0;
%     end
%     structFunctions(methodCounter) = constructFastSlowFunctions(['Varpro y = a*exp(b*x)+c*exp(d*x^f) ' axisPair], ...
%                 @(b,x) exp1(b(1:2),x)+expQ(b(3:5), x), B, exp1, B(1:2), expQ, B(3:5), 0, 0, r2);  
%     structFunctions(methodCounter).r2 = [r2, r2first, r2last];
% 
%     if graphState == 1
%         h = subplot(3,3,methodCounter);
%         annotateSubplot(h, 'Varpro y = a*exp(b*x)+c*exp(d*x^f)', gas, XAxisString, YAxisString, x, y, @(b,x)exp1(b(1:2),x)+expQ(b(3:5), x), B, exp1, B(1:2), expQ, B(3:5), r2);
%     end
%     if r2>1 || r2 <0 %If the data was not fit at all set random numbers for the next iteration
%         x0(:,7)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%     else
%         try
%             x0(:,7) = real(B);
%         catch
%             x0(:,7)    = [rand()*2, rand()*-1, rand()*2, rand()*-1, rand()];
%         end
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

function [Phi,dPhi,Ind] = adaex_exp(alpha,t)

% function [Phi,dPhi,Ind] = adaex(alpha,t)
% This is a sample user-defined function to be used by varpro.m.
%
% The model for this sample problem is
%
%   eta(t) = c1 exp(alpha1 t) + c2 exp(alpha2 t)
%
% Given t and alpha, we evaluate Phi, dPhi, and Ind.
%
% Dianne P. O'Leary and Bert W. Rust, September 2010.


% Evaluate Phi1 = exp(alpha1 t),
%          Phi2 = exp(alpha2 t),
% at each of the data points in t.

    Phi(:,1) = exp(alpha(1)*t);
    Phi(:,2) = exp(alpha(2)*t);

% The nonzero partial derivatives of Phi with respect to alpha are
%              d Phi_1 / d alpha_1 ,
%              d Phi_2 / d alpha_2 ,
% and this determines Ind.
% The ordering of the columns of Ind is arbitrary but must match dPhi.

    Ind = [1 2
           1 2];

% Evaluate the two nonzero partial derivatives of Phi at each of 
% the data points and store them in dPhi.

    dPhi(:,1) = t .* Phi(:,1);
    dPhi(:,2) = t .* Phi(:,2);
end  


function [Phi,dPhi,Ind] = adaex_pot(alpha,t)

    Phi(:,1) = exp(alpha(1)*t);
    Phi(:,2) = exp(alpha(2).^alpha(3)*t);

    Ind = [1 2 2
           1 2 3];

    dPhi(:,1) = t .* Phi(:,1);
    dPhi(:,2) = t .* alpha(3) .* alpha(2) .^ real(alpha(3) -1)   .* Phi(:,2);
    dPhi(:,3) = t .* alpha(2) .^ alpha(3) .* real(log(alpha(2))) .* Phi(:,2);
    
end
