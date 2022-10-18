function [] = time_stepper( pde_system, t, dt )

    switch pde_system.opts.timestep_method
        
        case 'FE'
            
            ForwardEuler( pde_system, t, dt );
        
        case 'SSPRK2'
            
            SSP_RK2( pde_system, t, dt );
            
        case 'SSPRK3'
            
            SSP_RK3( pde_system, t, dt );
            
        case {'BE','BE_LB'}
            
            BackwardEuler( pde_system, t, dt );
            
        case 'CN'
            
            CrankNicolson( pde_system, t, dt );

        case 'mM_FE'

            mM_ForwardEuler( pde_system, t, dt );

        case 'mM_SSPRK2'

            mM_SSPRK2( pde_system, t, dt )

        case 'mM_Dummy'

            mM_Dummy( pde_system, t, dt )
            
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

function [] = mM_ForwardEuler( pde_system, t, dt )

    i_rho_0 = 1;
    i_rho_1 = 2;
    i_rho_2 = 3;
    i_g     = 4;

    mM = MICRO_MACRO_1X1V( pde_system.opts, pde_system.unknowns{i_g}.dimensions );

    sv = pde_system.solution_vector;

    % --- Macro Moments at t^n ---

    rho_0 = {sv.unknowns{i_rho_0}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_0):sv.ubounds(i_rho_0)) ),...
             sv.unknowns{i_rho_1}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_1):sv.ubounds(i_rho_1)) ),...
             sv.unknowns{i_rho_2}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_2):sv.ubounds(i_rho_2)) )};

    % --- Maxwellian at t^n ---

    M_0 = sv.unknowns{i_g}.convert_to_wavelet( mM.evaluate_rhs_Maxwellian( pde_system.opts, rho_0, t ) );

    u = sv.copy_evolution_unknowns();

    % --- Evaluate RHS ---

    RHS = ComputeRHS_Explicit( pde_system, u, t );

    u = u + dt * RHS;

    sv.insert_evolution_unknowns( u );

    % --- Macro Moments at t^n+1 ---

    rho_1 = {sv.unknowns{i_rho_0}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_0):sv.ubounds(i_rho_0)) ),...
             sv.unknowns{i_rho_1}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_1):sv.ubounds(i_rho_1)) ),...
             sv.unknowns{i_rho_2}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_2):sv.ubounds(i_rho_2)) )};

    % --- Maxwellian at t^n+1 ---
    
    M_1 = sv.unknowns{i_g}.convert_to_wavelet( mM.evaluate_rhs_Maxwellian( pde_system.opts, rho_1, t + dt ) );
    
    % --- Add M^n - M^n+1 to g^n+1 ---

    sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        = sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        + pde_system.equations{i_g}.MultiplyInverseMassMatrix( pde_system.opts, M_0 - M_1, t + dt );

    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [] = mM_SSPRK2( pde_system, t, dt )

    i_rho_0 = 1;
    i_rho_1 = 2;
    i_rho_2 = 3;
    i_g     = 4;

    mM = MICRO_MACRO_1X1V( pde_system.opts, pde_system.unknowns{i_g}.dimensions );

    sv = pde_system.solution_vector;

    rho = {sv.unknowns{i_rho_0}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_0):sv.ubounds(i_rho_0)) ),...
           sv.unknowns{i_rho_1}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_1):sv.ubounds(i_rho_1)) ),...
           sv.unknowns{i_rho_2}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_2):sv.ubounds(i_rho_2)) )};

    M_n = sv.unknowns{i_g}.convert_to_wavelet( mM.evaluate_rhs_Maxwellian( pde_system.opts, rho, t ) );

    u_n = sv.copy_evolution_unknowns();

    RHS_n = ComputeRHS_Explicit( pde_system, u_n, t );

    u_1 = u_n + dt * RHS_n;

    sv.insert_evolution_unknowns( u_1 );

    rho = {sv.unknowns{i_rho_0}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_0):sv.ubounds(i_rho_0)) ),...
           sv.unknowns{i_rho_1}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_1):sv.ubounds(i_rho_1)) ),...
           sv.unknowns{i_rho_2}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_2):sv.ubounds(i_rho_2)) )};

    M_1 = sv.unknowns{i_g}.convert_to_wavelet( mM.evaluate_rhs_Maxwellian( pde_system.opts, rho, t ) );

    sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        = sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        - pde_system.equations{i_g}.MultiplyInverseMassMatrix( pde_system.opts, ( M_1 - M_n ), t + dt );

    u_1 = sv.copy_evolution_unknowns();

    RHS_1 = ComputeRHS_Explicit( pde_system, u_1, t+dt );

    u = u_n + 0.5 * dt * ( RHS_n + RHS_1 );

    sv.insert_evolution_unknowns( u );

    rho = {sv.unknowns{i_rho_0}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_0):sv.ubounds(i_rho_0)) ),...
           sv.unknowns{i_rho_1}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_1):sv.ubounds(i_rho_1)) ),...
           sv.unknowns{i_rho_2}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_2):sv.ubounds(i_rho_2)) )};

    M = sv.unknowns{i_g}.convert_to_wavelet( mM.evaluate_rhs_Maxwellian( pde_system.opts, rho, t ) );

    sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        = sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        - pde_system.equations{i_g}.MultiplyInverseMassMatrix( pde_system.opts, ( M - M_n ), t + dt );

    [ sv ] = EvaluateClosure( pde_system, sv, t + dt );

end

function [] = mM_Dummy( pde_system, t, dt )

    i_rho_0 = 1;
    i_rho_1 = 2;
    i_rho_2 = 3;
    i_g     = 4;

    mM = MICRO_MACRO_1X1V( pde_system.opts, pde_system.unknowns{i_g}.dimensions );

    sv = pde_system.solution_vector;

    rho = {sv.unknowns{i_rho_0}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_0):sv.ubounds(i_rho_0)) ),...
           sv.unknowns{i_rho_1}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_1):sv.ubounds(i_rho_1)) ),...
           sv.unknowns{i_rho_2}.convert_to_realspace( sv.fvec(sv.lbounds(i_rho_2):sv.ubounds(i_rho_2)) )};

    M = sv.unknowns{i_g}.convert_to_wavelet( mM.evaluate_rhs_Maxwellian( pde_system.opts, rho, t ) );

    u = sv.copy_evolution_unknowns();

    RHS = ComputeRHS_Explicit( pde_system, u, t );

    u = u + dt * RHS;

    sv.insert_evolution_unknowns( u );

    sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        = sv.fvec(sv.lbounds(i_g):sv.ubounds(i_g))...
        + pde_system.equations{i_g}.MultiplyInverseMassMatrix( pde_system.opts, M, t + dt );

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

    if strcmp( pde_system.opts.timestep_method, 'BE_LB' )

        Ax = @(x) BE_LHS_LB( pde_system, sv, x, t, dt );

    else

        Ax = @(x) BE_LHS( pde_system, sv, x, t, dt );

    end
    
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

    x_0 = sv.copy();
    
    [ sv.fvec, ~, relres, iter ] = bicgstabl( Ax, b, 1e-10, numel(b), [], [], x_0 );

    if ~pde_system.opts.quiet; disp(['    Backward Euler Solve (iter=',num2str(iter),')']); end

    assert( relres < 1e-7, 'BackwardEuler: bicgstabl failed' )

end

function [ Ax ] = BE_LHS( pde_system, sv, x, t, dt )

    opts = pde_system.opts;
    
    Ax = zeros(size(x));
    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        lo = sv.lbounds(i);
        hi = sv.ubounds(i);
        
        Ax(lo:hi) = equation.LHS_term.driver( opts, {x(lo:hi)}, t+dt );
        
        alpha = dt;
        if( strcmp( equation.type, 'closure' ) )
            alpha = 1.0;
        end
        
        for j = 1 : numel(equation.terms)
            
            term = equation.terms{j};
            
            Q = sv.get_input_unknowns( term.input_unknowns, x );
            
            Ax(lo:hi) = Ax(lo:hi) - alpha * term.driver( opts, Q, t+dt );
            
        end
        
    end
    
end

function [ Ax ] = BE_LHS_LB( pde_system, sv, x, t, dt )

    opts = pde_system.opts;
    
    Ax = zeros(size(x));
    for i = 1 : pde_system.num_eqs
        
        equation = pde_system.equations{i};
        
        lo = sv.lbounds(i);
        hi = sv.ubounds(i);
        
        Ax(lo:hi) = equation.LHS_term.driver( opts, {x(lo:hi)}, t+dt );
        
        alpha = dt;
        if( strcmp( equation.type, 'closure' ) )
            alpha = 1.0;
        end
        
        for j = 1 : numel(equation.terms)
            
            term = equation.terms{j};
            
            if     isequal( functions(term.descriptor{1}).function, '@(opts,Q,t)LB.evaluate_rhs_collision_operator_LB(opts,Q,t)' )
                
                num_Q = numel(term.input_unknowns);

                assert( num_Q == 5 )

                Q = cell(num_Q,1);

                Q{1} = sv.get_input_unknown( term.input_unknowns, 1, x );
                Q{2} = sv.get_input_unknown( term.input_unknowns, 2, x );
                Q{3} = sv.get_input_unknown( term.input_unknowns, 3 );
                Q{4} = sv.get_input_unknown( term.input_unknowns, 4 );
                Q{5} = sv.get_input_unknown( term.input_unknowns, 5 );

            elseif isequal( functions(term.descriptor{1}).function, '@(opts,Q,t)LB.evaluate_rhs_diffusion_LB(opts,Q,t)' )
                
                num_Q = numel(term.input_unknowns);

                assert( num_Q == 4 )

                Q = cell(num_Q,1);

                Q{1} = sv.get_input_unknown( term.input_unknowns, 1, x );
                Q{2} = sv.get_input_unknown( term.input_unknowns, 2 );
                Q{3} = sv.get_input_unknown( term.input_unknowns, 3 );
                Q{4} = sv.get_input_unknown( term.input_unknowns, 4 );

            else

                Q = sv.get_input_unknowns( term.input_unknowns, x );

            end
            
            Ax(lo:hi) = Ax(lo:hi) - alpha * term.driver( opts, Q, t+dt );
            
        end
        
    end
    
end