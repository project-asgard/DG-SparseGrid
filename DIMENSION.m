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
        
        function dim = DIMENSION(min,max,lev)
            
            if exist( 'min', 'var' )
                dim.min = min;
            end
            
            if exist( 'max', 'var' )
                dim.max = max;
            end
            
            if exist( 'lev', 'var' )
                dim.lev = lev;
            end
            
        end
        
    end
    
end
