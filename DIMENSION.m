classdef DIMENSION
    
    properties
        
        name = 'x';
        min = 0;
        max = 1;
        lev = 3;
        init_cond_fn = @(x,p) x.*0;
        jacobian = @(x,p,t) x.*0 + 1;
        
    end
    
    methods
        
        function dim = DIMENSION(min,max)
            
            if nargin == 2
               dim.min = min;
               dim.max = max;
            end
            
        end
        
    end
    
end
