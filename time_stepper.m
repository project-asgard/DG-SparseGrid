function [] = time_stepper( pde_system, t, dt )

    switch pde_system.opts.timestep_method
        
        case 'FE'
            
            ForwardEuler_NEW( pde_system, t, dt );
        
        case 'SSPRK2'
            
            SSP_RK2( pde_system, t, dt );
            
        case 'BE'
            
            BackwardEuler( pde_system, t, dt );
            
    end

end

function [] = ForwardEuler( pde_system, t, dt )

    num_eqs = pde_system.num_eqs;
    sv      = pde_system.solution_vector;
    
    RHS = zeros(size(sv.fvec));
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        equation.update_terms( pde_system.opts, t )
        
        if( strcmp(equation.type,'evolution') )
            
            lo = equation.unknown.lo_global;
            hi = equation.unknown.hi_global;
            
            for j = 1 : numel( equation.terms )
                
                RHS(lo:hi) = RHS(lo:hi) + equation.terms{j}.driver( pde_system.opts, sv, t );
                
            end
            
        end
        
    end
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp(equation.type,'evolution') )
            
            lo = equation.unknown.lo_global;
            hi = equation.unknown.hi_global;
            
            equation.InvertMassMatrix...
                ( pde_system.opts, sv, equation.LHS_term.driver( pde_system.opts, sv, t ) + dt * RHS(lo:hi), t );
            
        end
        
    end
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp( equation.type, 'closure' ) )
            
            equation.EvaluateClosure( pde_system.opts, sv, t );
            
        end
        
    end

end

function [] = ForwardEuler_NEW( pde_system, t, dt )

    sv = pde_system.solution_vector;

    u = sv.copy_evolution_unknowns();

    RHS = ComputeRHS( pde_system, u, t );
    
    u = u + dt * RHS;
    
    sv.insert_evolution_unknowns( u );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [] = SSP_RK2( pde_system, t, dt )

    sv = pde_system.solution_vector;

    u_0 = sv.copy_evolution_unknowns();

    RHS_0 = ComputeRHS( pde_system, u_0, t );
    
    u_1 = u_0 + dt * RHS_0;
    
    RHS_1 = ComputeRHS( pde_system, u_1, t );
    
    u_1 = u_0 + 0.5 * dt * ( RHS_0 + RHS_1 );
    
    sv.insert_evolution_unknowns( u_1 );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [ RHS ] = ComputeRHS( pde_system, u, t )

    sv_util = GLOBAL_SOLUTION_VECTOR_UTILITIES( );
    
    sv = sv_util.checkout_evolution( pde_system.solution_vector, u );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t );
    
    RHS = zeros(sv.size_evolution_unknowns(),1);
    
    os = 0;
    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp( equation.type, 'evolution' ) )
            
            equation.update_terms( pde_system.opts, t )
            
            unknown_size = equation.unknown.size();
            
            for j = 1 : numel( equation.terms )
                
                RHS(os+1:os+unknown_size)...
                    = RHS(os+1:os+unknown_size)...
                        + equation.terms{j}.driver( pde_system.opts, sv, t );
                
            end
            
            RHS(os+1:os+unknown_size)...
                = equation.MultiplyInverseMassMatrix( pde_system.opts, sv, RHS(os+1:os+unknown_size), t );
            
            os = os + unknown_size;
            
        end
        
    end
    
    sv.delete;

end

function [ sv ] = EvaluateClosure( pde_system, sv, t )

    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp( equation.type, 'closure' ) )
            
            equation.update_terms( pde_system.opts, t )
            
            equation.EvaluateClosure( pde_system.opts, sv, t );
            
        end
        
    end

end

function [] = BackwardEuler( pde_system, t, dt )

    num_eqs = pde_system.num_eqs;
    
    lo = pde_system.solution_vector.lbounds;
    hi = pde_system.solution_vector.ubounds;
    
    Ax = @(x) BE_LHS( pde_system, t, dt, x, lo, hi );
    
    b = zeros(size(pde_system.solution_vector.fvec));
    
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        equation.update_terms( pde_system.opts, t )
        
        if( strcmp( equation.type, 'evolution' ) )
            
            b(lo(i):hi(i)) = equation.LHS_term.driver( pde_system.opts, pde_system.solution_vector, t );
            
        end
        
    end
    
    [ x, ~, relres, iter ] = bicgstabl( Ax, b, 1e-10, numel(b) );
    
    assert( relres < 1e-7, 'BackwardEuler: bicgstabl failed' )
    
    for i = 1 : num_eqs
        
        pde_system.solution_vector.fvec(lo(i):hi(i)) = x(lo(i):hi(i));
        
    end

end

function [ Ax ] = BE_LHS( pde_system, t, dt, x, lo, hi )

    opts = pde_system.opts;
    sv   = pde_system.solution_vector;
    
    num_eqs = pde_system.num_eqs;
    
    % --- Hack ---
    tmp=sv.fvec;
    sv.fvec = x;
    % --- End Hack ---
    
    Ax = zeros(size(x));
    for i = 1 : num_eqs
        
        equation = pde_system.equations{i};
        
        Ax(lo(i):hi(i)) = equation.LHS_term.driver( opts, sv, t+dt, x(lo(i):hi(i)) );
        
        alpha = dt;
        if( strcmp(equation.type,'closure') )
            alpha = 1.0;
        end
        
        for j = 1 : numel(equation.terms)
            
            term = equation.terms{j};
            
            Ax(lo(i):hi(i)) = Ax(lo(i):hi(i)) - alpha * term.driver( opts, sv, t+dt );
            
        end
        
    end
    
    % --- Hack ---
    sv.fvec=tmp;
    % --- End Hack ---
    
end