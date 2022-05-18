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
    
    sv = pde_system.solution_vector;
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        equation.update_terms( pde_system.opts, t )
        
        if( strcmp(equation.type,'evolution') )
            
            RHS{i} = zeros(equation.unknown.size(),1);
            
            for j = 1 : numel( equation.terms )
                
                RHS{i} = RHS{i} + equation.terms{j}.driver( pde_system.opts, sv, t );
                
            end
            
        end
        
    end
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp(equation.type,'evolution') )
            
            equation.InvertMassMatrix...
                ( pde_system.opts, sv, equation.LHS_term.driver( pde_system.opts, sv, t ) + dt * RHS{i}, t );
            
        end
        
    end
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp(equation.type,'closure') )
            
            equation.EvaluateClosure( pde_system.opts, sv, t );
            
        end
        
    end

end

function [ RHS ] = ComputeRHS( pde_system, t )

    RHS = zeros(pde_system.solution_vector.size_evolution_unknowns(),1);

end

function [] = EvaluateClosure( )

end

function [] = BackwardEuler( pde_system, t, dt )

    num_eqs = pde_system.num_eqs;
    
    lo = uint64(zeros(num_eqs,1));
    hi = uint64(zeros(num_eqs,1));
    
    os = 0;
    for i = 1 : num_eqs
        lo(i) = os + 1;
        hi(i) = os + numel(pde_system.unknowns{i}.fval);
        os = hi(i);
    end
    
    Ax = @(x) BE_LHS( pde_system, t, dt, x, lo, hi );
    
    b = zeros(hi(end),1);
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        equation.update_terms( pde_system.opts, t )
        
        if( strcmp( equation.type, 'evolution' ) )
            
            b(lo(i):hi(i)) = equation.LHS_term.driver( pde_system.opts, t );
            
        end
        
    end
    
    [ x, ~, relres, iter ] = bicgstabl( Ax, b, 1e-10, numel(b) );
    
    assert( relres < 1e-7, 'BackwardEuler: bicgstabl failed' )
    
    for i = 1 : num_eqs
        
        pde_system.unknowns{i}.fval = x(lo(i):hi(i));
        
    end

end

function [ Ax ] = BE_LHS( pde_system, t, dt, x, lo, hi )

    opts = pde_system.opts;
    
    num_eqs = pde_system.num_eqs;
    
    x_cell = cell(num_eqs,1);
    for i = 1 : num_eqs
        
        x_cell{i} = x(lo(i):hi(i));
        
    end
    
    Ax = zeros(size(x));
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        Ax(lo(i):hi(i)) = equation.LHS_term.driver( opts, t+dt, x_cell(i) );
        
        alpha = dt;
        if( strcmp(equation.type,'closure') )
            alpha = 1.0;
        end
        
        for j = 1 : numel(equation.terms)
            
            term = equation.terms{j};
            
            Ax(lo(i):hi(i)) = Ax(lo(i):hi(i)) - alpha * term.driver( opts, t+dt, x_cell(term.input_g2l) );
            
        end
        
    end

end