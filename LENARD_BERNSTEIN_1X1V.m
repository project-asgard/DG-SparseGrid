classdef LENARD_BERNSTEIN_1X1V
    
    properties
        deg;
        dim_x;
        dim_v;
        lev_x;
        lev_v;
        BCs_X = 'PERIODIC';
        BCs_V = 'ZEROFLUX';
        T_min = 0.0;
        n_quad;
        alpha = 0.5;
        quad_x;
        quad_w;
        Phi;
        dPhi;
        PhiDn;
        PhiUp;
    end

    methods

        function obj = LENARD_BERNSTEIN_1X1V( opts, dimensions, BCs_X, BCs_V, T_min, n_quad, alpha_diffusion_flux )
            
            assert( numel(dimensions) == 2, 'Only Two Dimensions (1X and 1V)' )

            obj.deg   = opts.deg;
            obj.dim_x = dimensions{1};
            obj.lev_x = dimensions{1}.lev;
            obj.dim_v = dimensions{2};
            obj.lev_v = dimensions{2}.lev;

            if exist( 'BCs_X', 'var' ) && ~ isempty( BCs_X )
                
                obj.BCs_X = BCs_X;
                
            end
            
            if exist( 'BCs_V', 'var' ) && ~ isempty( BCs_V )
                
                obj.BCs_V = BCs_V;
                
            end

            if exist( 'T_min', 'var' ) && ~isempty( T_min )

                obj.T_min = T_min;

            end

            if exist( 'n_quad', 'var' ) && ~ isempty( n_quad )

                obj.n_quad = n_quad;

            else
            
                obj.n_quad = max( 10, 2 * obj.deg + 1 );

            end

            if exist( 'alpha_diffusion_flux', 'var' ) && ~isempty( alpha_diffusion_flux )

                obj.alpha = alpha_diffusion_flux;

            end

            [ obj.quad_x, obj.quad_w ] = lgwt( obj.n_quad, -1, +1 );
            
            obj.Phi   = lin_legendre ( obj.quad_x, obj.deg );
            obj.dPhi  = lin_dlegendre( obj.quad_x, obj.deg );
            obj.PhiDn = lin_legendre ( - 1.0     , obj.deg );
            obj.PhiUp = lin_legendre ( + 1.0     , obj.deg );

        end

        function rhs_diffusion = evaluate_rhs_diffusion_LB( obj, opts, Q, t )

            N_x   = 2^obj.lev_x;
            N_v   = 2^obj.lev_v;
            dof_x = obj.deg * N_x;
            dof_v = obj.deg * N_v;
            dof   = dof_x * dof_v;
            dof_K = obj.deg^2;

            dx = ( obj.dim_x.max - obj.dim_x.min ) / N_x;
            dv = ( obj.dim_v.max - obj.dim_v.min ) / N_v;

            i_f     = 1;
            i_rho_0 = 2;
            i_rho_1 = 3;
            i_rho_2 = 4;

            rho_0_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{i_rho_0}, [ obj.deg, N_x ] );
            rho_1_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{i_rho_1}, [ obj.deg, N_x ] );
            rho_2_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{i_rho_2}, [ obj.deg, N_x ] );

            [ ~, ~, T_q ] = Primitive( rho_0_q, rho_1_q, rho_2_q, obj.T_min );

            rhs_diffusion = zeros( dof, 1 );

            %% 2D quadrature Weights
            
            w_q = kron( obj.quad_w, obj.quad_w );
            
            %% 2D Basis Functions Evaluated in Quadrature Points on Local Element

            Phi_q = kron( obj.Phi, obj.Phi ) / sqrt( dx * dv );

            %% Defivative of 2D Basis Functions Evaluated in Quadrature Points on Local Element
            
            dPhi_q = kron( obj.Phi, obj.dPhi ) / sqrt( dx * dv ) / ( dv / 2 );

            % --- Volume Term ---
            
            for i_x = 1 : N_x
            for i_v = 1 : N_v

                i_K = ( i_x - 1 ) * N_v + i_v;
                
                f_q = Phi_q * Q{i_f}((i_K-1)*dof_K+1:i_K*dof_K);

                sqrtT_q = kron( sqrt( T_q(:,i_x) ), ones(obj.n_quad,1) );

                rhs_diffusion((i_K-1)*dof_K+1:i_K*dof_K)...
                    = - 0.25 * dx * dv * ( w_q .* sqrtT_q .* f_q )' * dPhi_q;

            end
            end
            
            % --- Numerical Flux ---

            PhiDn = kron( obj.Phi, obj.PhiDn ) / sqrt( dx * dv );
            PhiUp = kron( obj.Phi, obj.PhiUp ) / sqrt( dx * dv );

            F_num = zeros(obj.n_quad,N_v+1,N_x);

            for i_x = 1 : N_x
            for i_v = 1 : N_v + 1

                [ left  ] = GetIndicesV_L( i_x, N_x, i_v, N_v, dof_K, obj.BCs_V );
                [ right ] = GetIndicesV_R( i_x, N_x, i_v, N_v, dof_K, obj.BCs_V );

                f_L = PhiUp * Q{i_f}(left);
                f_R = PhiDn * Q{i_f}(right);

                F_num(:,i_v,i_x) = ( 1 - obj.alpha ) * f_L + obj.alpha * f_R;

            end
            end

            if( strcmp( obj.BCs_V, 'ZEROFLUX' ) )
                for i_x = 1 : N_x
                    F_num(:,1    ,i_x) = 0.0;
                    F_num(:,N_v+1,i_x) = 0.0;
                end
            end

            % --- Surface Term ---

            for i_x = 1 : N_x
            for i_v = 1 : N_v

                i_K = ( i_x - 1 ) * N_v + i_v;

                sqrtT_q = sqrt( T_q(:,i_x) );

                rhs_diffusion((i_K-1)*dof_K+1:i_K*dof_K)...
                    = rhs_diffusion((i_K-1)*dof_K+1:i_K*dof_K)...
                    + 0.5 * dx * (   PhiUp' * ( obj.quad_w .* sqrtT_q .* F_num(:,i_v+1,i_x) )...
                                   - PhiDn' * ( obj.quad_w .* sqrtT_q .* F_num(:,i_v  ,i_x) ) );

            end
            end

        end

        function rhs_collision_operator = evaluate_rhs_collision_operator_LB( obj, opts, Q, t )

            N_x   = 2^obj.lev_x;
            N_v   = 2^obj.lev_v;
            dof_x = obj.deg * N_x;
            dof_v = obj.deg * N_v;
            dof   = dof_x * dof_v;
            dof_K = obj.deg^2;

            dx = ( obj.dim_x.max - obj.dim_x.min ) / N_x;
            dv = ( obj.dim_v.max - obj.dim_v.min ) / N_v;

            i_f     = 1;
            i_q     = 2;
            i_rho_0 = 3;
            i_rho_1 = 4;
            i_rho_2 = 5;

            rho_0_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{i_rho_0}, [ obj.deg, N_x ] );
            rho_1_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{i_rho_1}, [ obj.deg, N_x ] );
            rho_2_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{i_rho_2}, [ obj.deg, N_x ] );

            [ ~, U_q, T_q ] = Primitive( rho_0_q, rho_1_q, rho_2_q, obj.T_min );

            rhs_collision_operator = zeros( dof, 1 );

            %% 2D quadrature Weights
            
            w_q = kron( obj.quad_w, obj.quad_w );
            
            %% 2D Basis Functions Evaluated in Quadrature Points on Local Element

            Phi_q = kron( obj.Phi, obj.Phi ) / sqrt( dx * dv );

            %% Defivative of 2D Basis Functions Evaluated in Quadrature Points on Local Element
            
            dPhi_q = kron( obj.Phi, obj.dPhi ) / sqrt( dx * dv ) / ( dv / 2 );

            % --- Volume Term ---
            
            for i_x = 1 : N_x
            for i_v = 1 : N_v

                i_K = ( i_x - 1 ) * N_v + i_v;

                v_q = obj.dim_v.min + ((i_v-1)+.5*(1.+obj.quad_x)) * dv;
                
                f_q = Phi_q * Q{i_f}((i_K-1)*dof_K+1:i_K*dof_K);
                q_q = Phi_q * Q{i_q}((i_K-1)*dof_K+1:i_K*dof_K);

                W_q     = kron( U_q(:,i_x)        , ones(obj.n_quad,1) )...
                        - kron( ones(obj.n_quad,1), v_q );
                sqrtT_q = kron( sqrt( T_q(:,i_x) ), ones(obj.n_quad,1) );

                rhs_collision_operator((i_K-1)*dof_K+1:i_K*dof_K)...
                    = 0.25 * dx * dv * ( w_q .* ( W_q .* f_q - sqrtT_q .* q_q ) )' * dPhi_q;

            end
            end

            % --- Numerical Flux ---

            PhiDn = kron( obj.Phi, obj.PhiDn ) / sqrt( dx * dv );
            PhiUp = kron( obj.Phi, obj.PhiUp ) / sqrt( dx * dv );

            F_num_f = zeros(obj.n_quad,N_v+1,N_x);
            F_num_q = zeros(obj.n_quad,N_v+1,N_x);

            for i_x = 1 : N_x
            for i_v = 1 : N_v + 1

                [ left  ] = GetIndicesV_L( i_x, N_x, i_v, N_v, dof_K, obj.BCs_V );
                [ right ] = GetIndicesV_R( i_x, N_x, i_v, N_v, dof_K, obj.BCs_V );

                % --- Drift Term ---

                f_L = PhiUp * Q{i_f}(left);
                f_R = PhiDn * Q{i_f}(right);

                W_q = U_q(:,i_x) - (obj.dim_v.min + (i_v-1) * dv);

                F_num_f(:,i_v,i_x) = 0.5 .* ( W_q .* f_L + W_q .* f_R - abs(W_q) .* ( f_R - f_L ) );

                % --- Diffusion Term ---

                q_L = PhiUp * Q{i_q}(left);
                q_R = PhiDn * Q{i_q}(right);

                F_num_q(:,i_v,i_x) = obj.alpha * q_L + ( 1 - obj.alpha ) * q_R;

            end
            end

            if( strcmp( obj.BCs_V, 'ZEROFLUX' ) )
                for i_x = 1 : N_x
                    F_num_f(:,1    ,i_x) = 0.0;
                    F_num_f(:,N_v+1,i_x) = 0.0;
                    F_num_q(:,1    ,i_x) = 0.0;
                    F_num_q(:,N_v+1,i_x) = 0.0;
                end
            end

            % --- Surface Term ---

            for i_x = 1 : N_x
            for i_v = 1 : N_v

                i_K = ( i_x - 1 ) * N_v + i_v;

                sqrtT_q = sqrt( T_q(:,i_x) );

                rhs_collision_operator((i_K-1)*dof_K+1:i_K*dof_K)...
                    = rhs_collision_operator((i_K-1)*dof_K+1:i_K*dof_K)...
                    - 0.5 * dx * (   PhiUp' * ( obj.quad_w .* ( F_num_f(:,i_v+1,i_x) - sqrtT_q .* F_num_q(:,i_v+1,i_x) ) )...
                                   - PhiDn' * ( obj.quad_w .* ( F_num_f(:,i_v  ,i_x) - sqrtT_q .* F_num_q(:,i_v  ,i_x) ) ) );

            end
            end

        end

    end
end

function [ indices ] = GetIndicesV_L( i_x, N_x, i_v, N_v, dof_K, BCs )

    i_K1   = (i_x-1)*N_v + 1;
    i_KN_v = (i_x-1)*N_v + N_v;

    first = (i_K1  -1)*dof_K+1:i_K1  *dof_K;
    last  = (i_KN_v-1)*dof_K+1:i_KN_v*dof_K;

    i_KL = (i_x-1)*N_v+i_v-1;

    if( i_v > 1 )
        indices = (i_KL-1)*dof_K+1:i_KL*dof_K;
    else
        switch BCs
            case 'PERIODIC'
                indices = last;
            case 'HOMOGENEOUS'
                indices = first;
            otherwise
                indices = first;
        end
    end

end

function [ indices ] = GetIndicesV_R( i_x, N_x, i_v, N_v, dof_K, BCs )

    i_K1   = (i_x-1)*N_v + 1;
    i_KN_v = (i_x-1)*N_v + N_v;

    first = (i_K1  -1)*dof_K+1:i_K1  *dof_K;
    last  = (i_KN_v-1)*dof_K+1:i_KN_v*dof_K;

    i_KR = (i_x-1)*N_v+i_v;

    if( i_v < N_v + 1 )
        indices = (i_KR-1)*dof_K+1:i_KR*dof_K;
    else
        switch BCs
            case 'PERIODIC'
                indices = first;
            case 'HOMOGENEOUS'
                indices = last;
            otherwise
                indices = last;
        end
    end

end

function [ D, U, T ] = Primitive( rho_0, rho_1, rho_2, T_min )

  D = rho_0;
  U = rho_1 ./ rho_0;
  T = max( ( 2.0 * rho_2 - rho_1.^2 ./ rho_0 ) ./ ( rho_0 ), T_min );

end