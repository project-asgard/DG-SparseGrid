classdef EULER_1D
    
    properties
        dim;
        lev;
        deg;
        BCs = 'PERIODIC';
        RiemannSolver = 'LLF';
        dim_v = 1;
        T_min = 0;
        quad_w;
        p_val;
        dp_val;
        p_L;
        p_R;
    end
    
    methods
        
        function euler_1D = EULER_1D( opts, dimensions, BCs, RiemannSolver, dim_v, T_min )
            
            assert( numel(dimensions) == 1, 'Only One Spatial Dimension' )
            
            euler_1D.dim = dimensions{1};
            euler_1D.lev = dimensions{1}.lev;
            euler_1D.deg = opts.deg;
            
            if exist( 'BCs', 'var' )
                
                euler_1D.BCs = BCs;
                
            end
            
            assert( any(strcmp(euler_1D.BCs,{'PERIODIC','HOMOGENEOUS'})),...
                    'BCs Must Be One of: PERIODIC,HOMOGENEOUS' )
            
            if exist( 'RiemannSolver', 'var' )
                
                euler_1D.RiemannSolver = RiemannSolver;
                
            end
            
            if exist( 'dim_v', 'var' )
                
                euler_1D.dim_v = dim_v;
                
            end
            
            if exist( 'T_min', 'var' )
                
                euler_1D.T_min = T_min;
                
            end
            
            [ quad_x, euler_1D.quad_w ] = lgwt( 2 * euler_1D.deg - 1, -1, +1 );
            
            euler_1D.p_val  = lin_legendre ( quad_x, euler_1D.deg );
            euler_1D.dp_val = lin_dlegendre( quad_x, euler_1D.deg );
            euler_1D.p_L    = lin_legendre ( - 1.0 , euler_1D.deg );
            euler_1D.p_R    = lin_legendre ( + 1.0 , euler_1D.deg );
            
        end
        
        function rhs_1 = evaluate_rhs_1( obj, opts, Q, t )
            
            dim = obj.dim;
            lev = obj.lev;
            deg = obj.deg;
            N   = 2^lev;
            dof = deg * N;
            h   = ( dim.max - dim.min ) / N;
            Jac = h / 2;
            
            quad_w = obj.quad_w;
            p_val  = obj.p_val  / sqrt( h );
            dp_val = obj.dp_val / sqrt( h ) / Jac;
            p_L    = obj.p_L    / sqrt( h );
            p_R    = obj.p_R    / sqrt( h );
            
            rhs_1 = zeros( dof, 1 );
            
            % --- Volume Term ---
            
            for i = 1 : N
                
                current = (i-1)*deg+1:i*deg;
                
                u_1 = p_val * Q{1}(current);
                u_2 = p_val * Q{2}(current);
                u_3 = p_val * Q{3}(current);
                
                [ D, U, T ] = Primitive( u_1, u_2, u_3, obj.dim_v, obj.T_min );
                
                F_1 = Flux_1( D, U, T, obj.dim_v );
                
                rhs_1(current) = ( dp_val' * ( quad_w .* F_1 ) ) * Jac;
                
            end
            
            % --- Numerical Fluxes ---
            
            F_1_Num = zeros(N+1,1);
            
            first = 1:deg;
            last  = (N-1)*deg+1:N*deg;
            for i = 1 : N + 1
                
                left  = GetIndices_L( i, deg, N, obj.BCs, first, last );
                right = GetIndices_R( i, deg, N, obj.BCs, first, last );
                
                % --- Left State ---
                
                u_1_L = p_R * Q{1}(left);
                u_2_L = p_R * Q{2}(left);
                u_3_L = p_R * Q{3}(left);
                
                [ D_L, U_L, T_L ] = Primitive( u_1_L, u_2_L, u_3_L, obj.dim_v, obj.T_min );
                
                Lambda_L = MaxAbsEigenvalue( D_L, U_L, T_L, obj.dim_v );
                
                F_1_L = Flux_1( D_L, U_L, T_L, obj.dim_v );
                
                % --- Right State ---
                
                u_1_R = p_L * Q{1}(right);
                u_2_R = p_L * Q{2}(right);
                u_3_R = p_L * Q{3}(right);
                
                [ D_R, U_R, T_R ] = Primitive( u_1_R, u_2_R, u_3_R, obj.dim_v, obj.T_min );
                
                Lambda_R = MaxAbsEigenvalue( D_R, U_R, T_R, obj.dim_v );
                
                F_1_R = Flux_1( D_R, U_R, T_R, obj.dim_v );
                
                alpha = max( [ Lambda_L, Lambda_R ] );
                
                F_1_Num(i) = NumericalFlux_LLF( u_1_L, u_1_R, F_1_L, F_1_R, alpha );
                
            end
            
            % --- Surface Term ---
            
            for i = 1 : N
                
                current = (i-1)*deg+1:i*deg;
                
                rhs_1(current)...
                    = rhs_1(current)...
                    - ( p_R' .* F_1_Num(i+1) - p_L' .* F_1_Num(i) );
                
            end
            
        end
        
        function rhs_2 = evaluate_rhs_2( obj, opts, Q, t )
            
            dim = obj.dim;
            lev = obj.lev;
            deg = obj.deg;
            N   = 2^lev;
            dof = deg * N;
            h   = ( dim.max - dim.min ) / N;
            Jac = h / 2;
            
            quad_w = obj.quad_w;
            p_val  = obj.p_val  / sqrt( h );
            dp_val = obj.dp_val / sqrt( h ) / Jac;
            p_L    = obj.p_L    / sqrt( h );
            p_R    = obj.p_R    / sqrt( h );
            
            rhs_2 = zeros( dof, 1 );
            
            % --- Volume Term ---
            
            for i = 1 : N
                
                current = (i-1)*deg+1:i*deg;
                
                u_1 = p_val * Q{1}(current);
                u_2 = p_val * Q{2}(current);
                u_3 = p_val * Q{3}(current);
                
                [ D, U, T ] = Primitive( u_1, u_2, u_3, obj.dim_v, obj.T_min );
                
                F_2 = Flux_2( D, U, T, obj.dim_v );
                
                rhs_2(current) = ( dp_val' * ( quad_w .* F_2 ) ) * Jac;
                
            end
            
            % --- Numerical Fluxes ---
            
            F_2_Num = zeros(N+1,1);
            
            first = 1:deg;
            last  = (N-1)*deg+1:N*deg;
            for i = 1 : N + 1
                
                left  = GetIndices_L( i, deg, N, obj.BCs, first, last );
                right = GetIndices_R( i, deg, N, obj.BCs, first, last );
                
                % --- Left State ---
                
                u_1_L = p_R * Q{1}(left);
                u_2_L = p_R * Q{2}(left);
                u_3_L = p_R * Q{3}(left);
                
                [ D_L, U_L, T_L ] = Primitive( u_1_L, u_2_L, u_3_L, obj.dim_v, obj.T_min );
                
                Lambda_L = MaxAbsEigenvalue( D_L, U_L, T_L, obj.dim_v );
                
                F_2_L = Flux_2( D_L, U_L, T_L, obj.dim_v );
                
                % --- Right State ---
                
                u_1_R = p_L * Q{1}(right);
                u_2_R = p_L * Q{2}(right);
                u_3_R = p_L * Q{3}(right);
                
                [ D_R, U_R, T_R ] = Primitive( u_1_R, u_2_R, u_3_R, obj.dim_v, obj.T_min );
                
                Lambda_R = MaxAbsEigenvalue( D_R, U_R, T_R, obj.dim_v );
                
                F_2_R = Flux_2( D_R, U_R, T_R, obj.dim_v );
                
                alpha = max( [ Lambda_L, Lambda_R ] );
                
                F_2_Num(i) = NumericalFlux_LLF( u_2_L, u_2_R, F_2_L, F_2_R, alpha );
                
            end
            
            % --- Surface Term ---
            
            for i = 1 : N
                
                current = (i-1)*deg+1:i*deg;
                
                rhs_2(current)...
                    = rhs_2(current)...
                    - ( p_R' .* F_2_Num(i+1) - p_L' .* F_2_Num(i) );
                
            end
            
        end
        
        function rhs_3 = evaluate_rhs_3( obj, opts, Q, t )
            
            dim = obj.dim;
            lev = obj.lev;
            deg = obj.deg;
            N   = 2^lev;
            dof = deg * N;
            h   = ( dim.max - dim.min ) / N;
            Jac = h / 2;
            
            quad_w = obj.quad_w;
            p_val  = obj.p_val  / sqrt( h );
            dp_val = obj.dp_val / sqrt( h ) / Jac;
            p_L    = obj.p_L    / sqrt( h );
            p_R    = obj.p_R    / sqrt( h );
            
            rhs_3 = zeros( dof, 1 );
            
            % --- Volume Term ---
            
            for i = 1 : N
                
                current = (i-1)*deg+1:i*deg;
                
                u_1 = p_val * Q{1}(current);
                u_2 = p_val * Q{2}(current);
                u_3 = p_val * Q{3}(current);
                
                [ D, U, T ] = Primitive( u_1, u_2, u_3, obj.dim_v, obj.T_min );
                
                F_3 = Flux_3( D, U, T, obj.dim_v );
                
                rhs_3(current) = ( dp_val' * ( quad_w .* F_3 ) ) * Jac;
                
            end
            
            % --- Numerical Fluxes ---
            
            F_3_Num = zeros(N+1,1);
            
            first = 1:deg;
            last  = (N-1)*deg+1:N*deg;
            for i = 1 : N + 1
                
                left  = GetIndices_L( i, deg, N, obj.BCs, first, last );
                right = GetIndices_R( i, deg, N, obj.BCs, first, last );
                
                % --- Left State ---
                
                u_1_L = p_R * Q{1}(left);
                u_2_L = p_R * Q{2}(left);
                u_3_L = p_R * Q{3}(left);
                
                [ D_L, U_L, T_L ] = Primitive( u_1_L, u_2_L, u_3_L, obj.dim_v, obj.T_min );
                
                Lambda_L = MaxAbsEigenvalue( D_L, U_L, T_L, obj.dim_v );
                
                F_3_L = Flux_3( D_L, U_L, T_L, obj.dim_v );
                
                % --- Right State ---
                
                u_1_R = p_L * Q{1}(right);
                u_2_R = p_L * Q{2}(right);
                u_3_R = p_L * Q{3}(right);
                
                [ D_R, U_R, T_R ] = Primitive( u_1_R, u_2_R, u_3_R, obj.dim_v, obj.T_min );
                
                Lambda_R = MaxAbsEigenvalue( D_R, U_R, T_R, obj.dim_v );
                
                F_3_R = Flux_3( D_R, U_R, T_R, obj.dim_v );
                
                alpha = max( [ Lambda_L, Lambda_R ] );
                
                F_3_Num(i) = NumericalFlux_LLF( u_3_L, u_3_R, F_3_L, F_3_R, alpha );
                
            end
            
            % --- Surface Term ---
            
            for i = 1 : N
                
                current = (i-1)*deg+1:i*deg;
                
                rhs_3(current)...
                    = rhs_3(current)...
                    - ( p_R' .* F_3_Num(i+1) - p_L' .* F_3_Num(i) );
                
            end
            
        end
        
    end
    
end

function [ indices ] = GetIndices_L( i, deg, N, BCs, first, last )

    if i == 1
        
        switch BCs
            case 'PERIODIC'
                indices = last;
            case 'HOMOGENEOUS'
                indices = first;
            otherwise
                indices = last;
        end
        
    else
        
        indices = (i-2)*deg+1:(i-1)*deg;
        
    end

end

function [ indices ] = GetIndices_R( i, deg, N, BCs, first, last )

    if i == N + 1
        
        switch BCs
            case 'PERIODIC'
                indices = first;
            case 'HOMOGENEOUS'
                indices = last;
            otherwise
                indices = first;
        end
        
    else
        
        indices = (i-1)*deg+1:i*deg;
        
    end

end

function [ D, U, T ] = Primitive( u_1, u_2, u_3, d_v, T_min )

  D = u_1;
  U = u_2 ./ u_1;
  T = max( ( 2.0 * u_3 - u_2.^2 ./ u_1 ) ./ ( d_v * u_1 ), T_min );

end

function [ F_1 ] = Flux_1( D, U, T, d_v )

  F_1 = D .* U;

end

function [ F_2 ] = Flux_2( D, U, T, d_v )

  F_2 = D .* ( U.^2 + T );

end

function [ F_3 ] = Flux_3( D, U, T, d_v )

  F_3 = 0.5 .* D .* ( U.^2 + (d_v+2) * T ) .* U;

end

function [ MaxAbsEig ] = MaxAbsEigenvalue( D, U, T, d_v )

    Cs = sqrt( ( 2.0 + d_v ) * T / d_v );
    MaxAbsEig = max( [ abs( U + Cs ), abs( U - Cs ) ] );

end

function [ F ] = NumericalFlux_LLF( u_L, u_R, F_L, F_R, alpha )

    F = 0.5 * ( F_L + F_R - alpha * ( u_R - u_L ) );

end