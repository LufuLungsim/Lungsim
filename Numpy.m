classdef Numpy

    properties (Access = private)
        
    end
    
    methods
        
       function c = legint(obj, a, cnt)
            
            c = a; %rollaxis(a, 0);
            
            k = zeros(1, cnt);
            for i = 0 : cnt-1 
                n = size(c, 1);
                if n == 1 && sum(c(0+1,:) == 0) == length(c(0+1,:))
                    c(0+1, :) = c(0+1,:) + k(i+1);
                else
                    tmp = zeros(n+1, size(c, 2)); 
                    tmp(0+1, :) = c(0+1, :)*0;
                    tmp(1+1, :) = c(0+1, :);
                end
                
                if n > 1
                    tmp(2+1, :) = c(1+1, :)/3;
                end
                
                for j = 2 : n-1
                    t = c(j+1, :)/(2*j + 1);
                    tmp(j + 1+1,:) = t;
                    tmp(j - 1+1,:) = tmp(j - 1+1,:) - t;
                end
                tmp(0+1, :) = tmp(0+1, :) + k(1, i+1) - obj.legval(0, tmp);
                c = tmp;
            end
            %c = np.rollaxis(c, 0, 1);
        end

        function v = legval(obj, x, c)
           
            %c = np.array(c, ndmin=1, copy=0)
            
            if length(c) == 1
                c0 = c(0+1, :);
                c1 = 0;
            else
                if size(c, 1) == 2
                    c0 = c(0+1,:);
                    c1 = c(1+1,:);
                else
                    nd = size(c, 1);
                    c0 = c(end - 1,:);
                    c1 = c(end, :);
                    for i = 3 : size(c, 1)
                        tmp = c0;
                        nd = nd - 1;
                        c0 = c(end-(i) +1,:) - (c1*(nd - 1))/nd;
                        c1 = tmp + (c1.*x*(2*nd - 1))/nd;
                    end
                end
            end
            v = c0 + c1.*x;       
        end
        
    end
    
end
