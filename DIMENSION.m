classdef DIMENSION
    
    properties
        
        name = 'x';
        min = 0;
        max = 1;
        lev = 3;
        init_cond_fn = @(x,p) x.*0;
        moment_dV = @(x,p,t,d) x.*0 + 1;
        mass_mat = [];
        
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
