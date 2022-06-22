classdef PDE_SYSTEM < handle
    
    properties
        opts;
        num_eqs;
        unknowns;
        equations;
        solution_vector;
        solution_vector_utilities;
        set_dt = @(pde_system,CFL) 1;
    end
    
    methods
        
        function pde_system = PDE_SYSTEM( opts, equations, set_dt )
 
            pde_system.opts      = opts;
            pde_system.equations = equations;
            pde_system.set_dt    = set_dt;
            
            pde_system.num_eqs = numel( pde_system.equations );
            
            pde_system.unknowns = cell( pde_system.num_eqs, 1 );
            for i = 1 : pde_system.num_eqs
                
                pde_system.unknowns{i} = pde_system.equations{i}.unknown;
                
            end
            
            pde_system.solution_vector...
                = GLOBAL_SOLUTION_VECTOR( pde_system.unknowns );
            
            pde_system.solution_vector_utilities...
                = GLOBAL_SOLUTION_VECTOR_UTILITIES( );
            
        end
        
        function set_initial_conditions( obj, t )
            
            sv = obj.solution_vector;
            for i = 1 : numel( obj.unknowns )
                
                sv.fvec(sv.lbounds(i):sv.ubounds(i))...
                    = obj.unknowns{i}.get_initial_conditions( obj.opts, t );
                
            end
            
        end
        
        function [ error ] = compute_error( obj, t )
            
            num_unknowns = numel( obj.unknowns );
            
            error = zeros(num_unknowns,1);
            
            for i = 1 : num_unknowns
                
                sv = obj.solution_vector;
                
                error(i) = obj.unknowns{i}.compute_error( obj.opts, sv.fvec(sv.lbounds(i):sv.ubounds(i)), t );
                
            end
            
        end
        
    end
    
end

