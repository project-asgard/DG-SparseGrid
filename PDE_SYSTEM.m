classdef PDE_SYSTEM < handle
    
    properties
        opts;
        num_eqs;
        unknowns;
        equations;
        solution_vector;
        solution_vector_utilities;
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
            
            pde_system.solution_vector...
                = GLOBAL_SOLUTION_VECTOR( pde_system.unknowns );
            
            pde_system.solution_vector_utilities...
                = GLOBAL_SOLUTION_VECTOR_UTILITIES( );
            
            pde_system.create_global_to_local_map;
            
        end
        
        function set_initial_conditions( obj )
            
            sv = obj.solution_vector;
            for i = 1 : numel( obj.unknowns )
                
                sv.fvec(sv.lbounds(i):sv.ubounds(i))...
                    = obj.unknowns{i}.get_initial_conditions( obj.opts );
                
            end
            
        end
        
        function create_global_to_local_map( obj )
            
            num_unknowns = numel( obj.unknowns );
            for i = 1 : obj.num_eqs
                
                equation = obj.equations{i};
                
                for j = 1 : numel( equation.terms )
                    
                    term = equation.terms{j};
                    
                    num_input_unknowns = numel( term.input_unknowns );
                    
                    term.input_g2l = int32(zeros(size(num_input_unknowns,1)));
                    
                    for k = 1 : num_input_unknowns
                        
                        for l = 1 : num_unknowns
                            
                            if( obj.unknowns{l} == term.input_unknowns{k} )
                                
                                term.input_g2l(k) = l;
                                
                            end
                            
                        end
                        
                    end
                    
                end
                
            end
            
        end
        
    end
    
end

