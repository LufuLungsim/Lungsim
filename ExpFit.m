classdef ExpFit

    properties (Access = private)

        t_;
        z_;
        m_;
        p_; 
        
        Legendre_;
        trapezoidWeights_;       
        
        L_;
        odeSolve_;
        yhat_;
        x_;
    end
    
    methods
        
        function [lambda, coefficients, zFit] = do(obj, t, z, m, p)
            
            obj.t_ = t;
            obj.z_ = z;
            obj.m_ = m;
            obj.p_ = p;
            
            n = length(t);

            dt = t(2:end) - t(1:end-1);
            scale = 2 / (t(end) - t(1));

            trapezoidWeights = zeros(1, n);
            trapezoidWeights(1) = 0.5 * dt(1);
            trapezoidWeights(2:end-1) = 0.5 * (dt(2:end) + dt(1:end-1));
            trapezoidWeights(end) = 0.5 * dt(end);
            trapezoidWeights = trapezoidWeights * scale;
            obj.trapezoidWeights_ = trapezoidWeights;

            x = 2 * (t - t(1)) / (t(end) - t(1)) - 1;
            obj.x_ = x;

            Legendre = zeros(m, n);
            Legendre(1, :) = 1;
            Legendre(2, :) = x;

            for k = 2 : m - 1
                Legendre(k + 1, :) = ((2*k-1) * Legendre(k, :) .* x'  - (k - 1) * Legendre(k - 1, :)) / k;
            end
            
            obj.Legendre_ = Legendre;
            obj.L_ = Legendre(1:p, :)';
            
            
            odeSolve = OdeSolve;
            odeSolve = odeSolve.init(m, p);
            obj.odeSolve_ = odeSolve;
            
            a0 = ones(p+1, 1);
            options = solveLMF('default');
            options = solveLMF(options,'XTol',1e-6,'FTol',1e-12,'ScaleD',1/scale,'Display',0);
            residualFunction = @(a) obj.residual(a);
            warning('off','MATLAB:nearlySingularMatrix')
            warning('off','MATLAB:illConditionedMatrix')
            warning('off','MATLAB:singularMatrix')
            [aopt1]=solveLMF(residualFunction,a0,options);
            warning('on','MATLAB:nearlySingularMatrix')
            warning('on','MATLAB:illConditionedMatrix')
            warning('on','MATLAB:singularMatrix')

            aoptReverse = fliplr(aopt1')';
            lambda = roots(aoptReverse) * scale;
            
            warning('off','MATLAB:rankDeficientMatrix')
            coefficients = exp(obj.t_ * [ 0; lambda]') \ obj.z_;
            warning('on','MATLAB:rankDeficientMatrix')
            
            zFit = exp(obj.t_ * [ 0; lambda]') * coefficients;
        end
    end
    
    methods (Access = private)
        
        function yhat = getYhat(obj, a)
           
            obj.odeSolve_ = obj.odeSolve_.setCoeffs(a);
            v = zeros(obj.m_, 1);
            v(1) = 1;
            c = diag(((0:obj.m_-1) + 0.5)) *  obj.Legendre_ * (obj.trapezoidWeights_' .* obj.z_);
            obj.odeSolve_ = obj.odeSolve_.solve(c(1:obj.p_), v);
            
            yhat = obj.odeSolve_.getY()';
        end
        
        function res = residual(obj, a)
            
            obj.odeSolve_ = obj.odeSolve_.setCoeffs(a);
            v = zeros(obj.m_, 1);
            v(1) = 1;
            c = diag(((0:obj.m_-1) + 0.5)) *  obj.Legendre_ * (obj.trapezoidWeights_' .* obj.z_);
            obj.odeSolve_ = obj.odeSolve_.solve(c(1:obj.p_), v);
            
            obj.yhat_ = obj.odeSolve_.getY()';
            
            numpy = Numpy;
            y = numpy.legval(obj.x_, obj.yhat_);
            
            res = obj.z_ - y;
        end
    end
end

