classdef PDE_SYSTEM < handle
    
    properties
        opts;
        num_eqs;
        unknowns;
        equations;
    end
    
    methods
        
        function pde_system = PDE_SYSTEM( opts, equations )
 
            pde_system.opts      = opts;
            pde_system.equations = equations;
            
            pde_system.num_eqs = numel( pde_system.equations );
            
            unknowns = cell( pde_system.num_eqs );
            for i = pde_system.num_eqs
                
                unknowns{i} = pde_system.equations{i}.unknown;
                
            end
            
        end
        
    end
end

