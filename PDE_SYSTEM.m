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
        
        function lev_vec = get_lev_vec(obj)
            num_dims = numel(obj.pde_solutions{1}.dimensions); % This must be changed for multiple solutions!
            for d=1:num_dims
                lev_vec(d,1) = obj.pde_solutions{1}.dimensions{d}.lev; % This must be changed for multiple solutions!
            end
        end
        
        function initialize( obj )
            
            [obj.hash_table.elements, obj.hash_table.elements_idx]...
                = hash_table_sparse_nD( obj.get_lev_vec, obj.opts.max_lev, obj.opts.grid_type );
            
            for i = 1 : obj.num_comps
                
                obj.pde_solutions{i}.set_initial_conditions( obj.opts, obj.hash_table )
                
                obj.pde_solutions{i}.initialize_auxiliary()
                
            end
            
        end
        
    end
end

