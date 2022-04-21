classdef PDE_SYSTEM < handle
    
    properties
        opts;
        hash_table;
        num_comps = 0;
        pde_solutions = {};
        pde_terms = {};
    end
    
    methods
        
        function pde_system = PDE_SYSTEM( opts, pde_solutions, pde_terms )
 
            pde_system.opts          = OPTS( opts );
            pde_system.num_comps     = numel(pde_solutions);
            pde_system.pde_solutions = pde_solutions;
            pde_system.pde_terms     = pde_terms;
            
        end
        
        
        
        function initialize( obj )
            
            
            
            for i = 1 : obj.num_comps
                
                obj.pde_solutions{i}.set_initial_conditions( obj.opts, obj.hash_table )
                
                obj.pde_solutions{i}.initialize_auxiliary()
                
            end
            
        end
        
    end
end

