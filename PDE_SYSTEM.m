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
            
            pde_system.unknowns = cell( pde_system.num_eqs, 1 );
            for i = 1 : pde_system.num_eqs
                
                pde_system.unknowns{i} = pde_system.equations{i}.unknown;
                
            end
            
        end
        
        function [ fvals ] = get_fvals( obj )
            
            for i = 1 : numel(obj.unknowns)
                
                fvals{i} = obj.unknowns{i}.fval;
                
            end
            
        end
        
    end
end

