classdef OdeSolve

    properties (Access = private)
        
        m_;
        p_;
        B_;
        a_;
        c_;
        y_;

        A_;
        
    end
    
    methods
        
        function obj = init(obj, m, p)
            
            obj.m_ = m;
            obj.p_ = p;
            obj.a_ = zeros(1, p);
            obj.c_ = zeros(1, p);
            obj.y_ = zeros(1, m);
            obj.B_ = {};
            
            for k = 0 : p-1
                matLegendre = obj.intMatLegendre(m, p-k);
                matLegendreBlock = matLegendre(p+1:(end-(p-k)),:);              
                tempMatrix = zeros(1, 1);
                tempMatrix( 1 : size(matLegendreBlock, 1), 1 : size(matLegendreBlock, 2)) = matLegendreBlock;
                obj.B_{p+1 - k} = tempMatrix;                
            end
            
            matLegendre = obj.intMatLegendre(m, 0);
            matLegendreBlock = matLegendre(p+1:end,:);
            tempMatrix = zeros(1, 1);
            tempMatrix( 1 : size(matLegendreBlock, 1), 1 : size(matLegendreBlock, 2)) = matLegendreBlock;
            obj.B_{1} = tempMatrix;
            
        end
        
        function obj = setCoeffs(obj, a)
            
            obj.a_ = a;
            A = zeros(obj.m_ - obj.p_, obj.m_);

            for k = 0 : obj.p_
               A = A + obj.a_(k+1) * obj.B_{obj.p_ - k + 1};
            end
            
            obj.A_ = A;
            
        end
        
        function obj = solve(obj, c, f)
            
            obj.c_ = c;
            obj.y_(1:obj.p_) = obj.c_;
          
            cond = obj.A_(:, 1:obj.p_);
            
            rhs = obj.B_{end} * f - cond*c;
            
            temp = obj.A_(:,obj.p_+1:end) \ rhs;
            obj.y_(obj.p_+1:end) = temp;
           
        end
        
        function y = getY(obj)
           
            y = obj.y_;
            
        end
        
%         function tangent = getTangent(obj)
%            
%             dy = zeros(obj.m_, obj.p_+1);
%             
%             for k = 0 : obj.p_
%                By = obj.B_{obj.p_ - k + 1} * obj.y_';
%                temp = obj.A_ \ By;
%                dy(obj.p_+1:end, k+1) = -temp(obj.p_+ 1:end);          
%             end
%             
%             tangent = dy;
%             
%         end
        
    end
    
    methods (Access = private)
        
        function B = intMatLegendre(obj, n, order)
            numpy = Numpy;            
            I = eye(n);
            B = numpy.legint(I, order);
        end
        
    end
    
end
