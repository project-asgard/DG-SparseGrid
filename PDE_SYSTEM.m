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
            
            % --- Create global to local map ---
            
            num_unknowns = numel( pde_system.unknowns );
            for i = 1 : pde_system.num_eqs
                
                equation = pde_system.equations{i};
                
                for j = 1 : numel( equation.terms )
                    
                    term = equation.terms{j};
                    
                    num_input_unknowns = numel( term.input_unknowns );
                    
                    term.input_g2l = int32(zeros(size(num_input_unknowns,1)));
                    
                    for k = 1 : num_input_unknowns
                        
                        for l = 1 : num_unknowns
                            
                            if( pde_system.unknowns{l} == term.input_unknowns{k} )
                                
                                term.input_g2l(k) = l;
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
        function [ fvals ] = get_fvals( obj )
            
            for i = 1 : numel(obj.unknowns)
                
                fvals{i} = obj.unknowns{i}.fval;
                
            end
            
        end
        
    end
end

