function [] = time_stepper( pde_system, t, dt )

    switch pde_system.opts.timestep_method
        
        case 'FE'
            
            ForwardEuler( pde_system, t, dt );
        
        case 'SSPRK2'
            
            SSP_RK2( pde_system, t, dt );
            
        case 'SSPRK3'
            
            SSP_RK3( pde_system, t, dt );
            
        case 'BE'
            
            BackwardEuler( pde_system, t, dt );
            
        case 'CrankNicolson'
            
            CrankNicolson( pde_system, t, dt );
            
    end

end

function [] = ForwardEuler( pde_system, t, dt )

    sv = pde_system.solution_vector;

    u = sv.copy_evolution_unknowns();

    RHS = ComputeRHS_Explicit( pde_system, u, t );
    
    u = u + dt * RHS;
    
    sv.insert_evolution_unknowns( u );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [] = SSP_RK2( pde_system, t, dt )

    sv = pde_system.solution_vector;

    u_0 = sv.copy_evolution_unknowns();

    RHS_0 = ComputeRHS_Explicit( pde_system, u_0, t );
    
    u_1 = u_0 + dt * RHS_0;
    
    RHS_1 = ComputeRHS_Explicit( pde_system, u_1, t + dt );
    
    u_1 = u_0 + 0.5 * dt * ( RHS_0 + RHS_1 );
    
    sv.insert_evolution_unknowns( u_1 );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [] = SSP_RK3( pde_system, t, dt )

    sv = pde_system.solution_vector;

    u_0 = sv.copy_evolution_unknowns();

    RHS_0 = ComputeRHS_Explicit( pde_system, u_0, t );
    
    u_1 = u_0 + dt * RHS_0;
    
    RHS_1 = ComputeRHS_Explicit( pde_system, u_1, t + dt );
    
    u_1 = u_0 + 0.25 * dt * ( RHS_0 + RHS_1 );
    
    RHS_2 = ComputeRHS_Explicit( pde_system, u_1, t + 0.25 * dt );
    
    u_1 = u_0 + (2.0/3.0) * dt * ( 0.25 * ( RHS_0 + RHS_1 ) + RHS_2 );
    
    sv.insert_evolution_unknowns( u_1 );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [] = BackwardEuler( pde_system, t, dt )

    sv = pde_system.solution_vector;

    u = sv.copy_evolution_unknowns();

    RHS = ComputeRHS_Implicit( pde_system, u, t, dt );
    
    u = u + dt * RHS;
    
    sv.insert_evolution_unknowns( u );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [] = CrankNicolson( pde_system, t, dt )

    sv = pde_system.solution_vector;

    u = sv.copy_evolution_unknowns();

    RHS = ComputeRHS_Explicit( pde_system, u, t );
    
    u = u + 0.5 * dt * RHS;
    
    RHS = ComputeRHS_Implicit( pde_system, u, t + 0.5 * dt, 0.5 * dt );
    
    u = u + 0.5 * dt * RHS;
    
    sv.insert_evolution_unknowns( u );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [ RHS ] = ComputeRHS_Explicit( pde_system, u, t )

    sv_util = GLOBAL_SOLUTION_VECTOR_UTILITIES( );
    
    sv = sv_util.checkout_evolution( pde_system.solution_vector, u );
    
    [ sv ] = EvaluateClosure( pde_system, sv, t );
    
    RHS = EvaluateRHS( pde_system, sv, t );
    
    sv.delete;

end

function [ RHS ] = ComputeRHS_Implicit( pde_system, u, t, dt )

    if( dt == 0.0 )
        
        RHS = zeros( size( u ) );
        
    else
        
        RHS = - u;
        
        sv_util = GLOBAL_SOLUTION_VECTOR_UTILITIES( );
        
        sv = sv_util.checkout_evolution( pde_system.solution_vector, u );
    
        [ sv ] = EvaluateClosure( pde_system, sv, t );
        
        BackwardEuler_Solve( pde_system, sv, t, dt );
        
        RHS = ( RHS + sv.copy_evolution_unknowns() ) ./ dt;
        
        sv.delete;
        
    end

end

function [ sv ] = EvaluateClosure( pde_system, sv, t )

    % --- Assumes Closure Unknowns Only Depend on Evolution Unknowns ---

    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp( equation.type, 'closure' ) )
            
            equation.update_terms( pde_system.opts, t )
            
            lo = equation.unknown.lo_global;
            hi = equation.unknown.hi_global;
            
            sv.fvec(lo:hi) = 0.0;
            for j = 1 : numel( equation.terms )
                
                term = equation.terms{j};
                
                Q = sv.get_input_unknowns( term.input_unknowns );
                
                sv.fvec(lo:hi) = sv.fvec(lo:hi) + term.driver( pde_system.opts, Q, t );
                
            end
            
            sv.fvec(lo:hi)...
                = equation.MultiplyInverseMassMatrix( pde_system.opts, sv.fvec(lo:hi), t );
            
        end
        
    end

end

function [ RHS ] = EvaluateRHS( pde_system, sv, t )

    RHS = sv.zeros_evolution();
    
    os = 0;
    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        if( strcmp( equation.type, 'evolution' ) )
            
            equation.update_terms( pde_system.opts, t )
            
            unknown_size = equation.unknown.size();
            
            for j = 1 : numel( equation.terms )
                
                Q = sv.get_input_unknowns( equation.terms{j}.input_unknowns );
                
                RHS(os+1:os+unknown_size)...
                    = RHS(os+1:os+unknown_size)...
                        + equation.terms{j}.driver( pde_system.opts, Q, t );
                
            end
            
            RHS(os+1:os+unknown_size)...
                = equation.MultiplyInverseMassMatrix( pde_system.opts, RHS(os+1:os+unknown_size), t );
            
            os = os + unknown_size;
            
        end
        
    end

end

function [] = BackwardEuler_Solve( pde_system, sv, t, dt )

    Ax = @(x) BE_LHS( pde_system, x, t, dt );
    
    b = sv.zeros();
    
    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        equation.update_terms( pde_system.opts, t )
        
        if( strcmp( equation.type, 'evolution' ) )
            
            lo = sv.lbounds(i);
            hi = sv.ubounds(i);
            
            b(lo:hi) = equation.LHS_term.driver( pde_system.opts, {sv.fvec(lo:hi)}, t );
            
        end
        
    end
    
    [ sv.fvec, ~, relres, iter ] = bicgstabl( Ax, b, 1e-10, numel(b) );
    
    assert( relres < 1e-7, 'BackwardEuler: bicgstabl failed' )

end

function [ Ax ] = BE_LHS( pde_system, x, t, dt )

    opts = pde_system.opts;
    
    Ax = zeros(size(x));
    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        lo = pde_system.solution_vector.lbounds(i);
        hi = pde_system.solution_vector.ubounds(i);
        
        Ax(lo:hi) = equation.LHS_term.driver( opts, {x(lo:hi)}, t+dt );
        
        alpha = dt;
        if( strcmp( equation.type, 'closure' ) )
            alpha = 1.0;
        end
        
        for j = 1 : numel(equation.terms)
            
            term = equation.terms{j};
            
            Q = pde_system.solution_vector.get_input_unknowns( term.input_unknowns, x );
            
            Ax(lo:hi) = Ax(lo:hi) - alpha * term.driver( opts, Q, t+dt );
            
        end
        
    end
    
end