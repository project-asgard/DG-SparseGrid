function [] = time_stepper( pde_system, t, dt )

    switch pde_system.opts.timestep_method
        
        case 'FE'
            
            ForwardEuler( pde_system, t, dt );
            
        case 'BE'
            
            BackwardEuler( pde_system, t, dt );
            
    end

end

function [] = ForwardEuler( pde_system, t, dt )

    num_eqs = pde_system.num_eqs;
    
    RHS = cell( num_eqs, 1 );
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp(equation.type,'evolution') )
            
            RHS{i} = zeros(size(equation.unknown.fval));
            
            for j = 1 : numel( equation.terms )
                
                RHS{i} = RHS{i} + equation.terms{j}.driver( pde_system.opts, t );
                
            end
            
        end
        
    end
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp(equation.type,'evolution') )
            
            equation.InvertMassMatrix...
              ( pde_system.opts, t, equation.LHS_term.driver(pde_system.opts,t) + dt * RHS{i} );
            
        end
        
    end
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp(equation.type,'closure') )
            
            equation.EvaluateClosure( pde_system.opts, t );
            
        end
        
    end

end

function [] = BackwardEuler( pde_system, t, dt )

    

end

function [ Ax ] = BE_LHS( pde_system, t, dt, x )

    num_eqs = pde_system.num_eqs;
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        for j = 1 : numel(equation.terms)
            
            
            
        end
        
    end

end