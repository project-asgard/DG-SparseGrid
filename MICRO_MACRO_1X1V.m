classdef MICRO_MACRO_1X1V
    
    properties
        deg;
        dim_x;
        dim_v;
        lev_x;
        lev_v;
        BCs_X = 'PERIODIC';
        BCs_V = 'PERIODIC';
        T_min = 0.0;
        n_quad;
        quad_x;
        quad_w;
        Phi;
        dPhi;
        PhiDn;
        PhiUp;
    end
    
    methods
        
        function obj = MICRO_MACRO_1X1V( opts, dimensions, BCs_X, BCs_V, T_min, n_quad )
            
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
            
                obj.n_quad = max( 10, 2 * obj.deg + 1 ); % --- May need adjustment to accurately integrate Maxwellian in v.

            end
            
            [ obj.quad_x, obj.quad_w ] = lgwt( obj.n_quad, -1, +1 );
            
            obj.Phi   = lin_legendre ( obj.quad_x, obj.deg );
            obj.dPhi  = lin_dlegendre( obj.quad_x, obj.deg );
            obj.PhiDn = lin_legendre ( - 1.0     , obj.deg );
            obj.PhiUp = lin_legendre ( + 1.0     , obj.deg );
            
        end
        
        function rhs_Maxwellian = evaluate_rhs_Maxwellian( obj, opts, Q, t )
            
            tic
            
            %% Evaluates int_K M(rho_h) phi_h dvdx
            %  K = [ v_L, v_H ] x [ x_L, x_H ] (local phase-space element)
            %  M is Maxwellian evaluated with moments rho_h
            %  phi_h are test functions
            
            N_x   = 2^obj.lev_x;
            N_v   = 2^obj.lev_v;
            dof_x = obj.deg * N_x;
            dof_v = obj.deg * N_v;
            dof   = dof_x * dof_v;
            dof_K = obj.deg^2;
            
            dx = ( obj.dim_x.max - obj.dim_x.min ) / N_x;
            dv = ( obj.dim_v.max - obj.dim_v.min ) / N_v;
            
            rho_1_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{1}, [ obj.deg, N_x ] );
            rho_2_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{2}, [ obj.deg, N_x ] );
            rho_3_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{3}, [ obj.deg, N_x ] );
            
            [ D_q, U_q, T_q ] = Primitive( rho_1_q, rho_2_q, rho_3_q, obj.T_min );
            
            rhs_Maxwellian = zeros( dof, 1 );
            
            %% 2D quadrature Weights
            
            w_q = kron( obj.quad_w, obj.quad_w ); % Move to initialization?
            
            %% 2D Basis Functions Evaluated in Quadrature Points on Local Element

            Phi_q = kron( obj.Phi, obj.Phi ) / sqrt( dx * dv ); % Move to initialization?
            
            i_K = 0;
            for i_x = 1 : N_x
            for i_v = 1 : N_v
                
                i_K = i_K + 1;
                
                v_q = obj.dim_v.min + ((i_v-1)+.5*(1.+obj.quad_x)) * dv;
                
                M_q = Maxwellian_q( D_q(:,i_x), U_q(:,i_x), T_q(:,i_x), v_q );
                
                rhs_Maxwellian((i_K-1)*dof_K+1:i_K*dof_K)...
                    = 0.25 * dx * dv * ( ( w_q .* M_q )' * Phi_q );
                
            end
            end
            
            toc
            
        end
        
        function rhs_vDotGradMaxwellian = evaluate_rhs_vDotGradMaxwellian( obj, opts, Q, t )
            
            tic
            
            %% Evaluates int_K (vM(rho_h))_x phi_h dvdx
            %  K = [ v_L, v_H ] x [ x_L, x_H ] (local phase-space element)
            %  M is Maxwellian evaluated with moments rho_h
            %  phi_h are test functions
            
            N_x   = 2^obj.lev_x;
            N_v   = 2^obj.lev_v;
            dof_x = obj.deg * N_x;
            dof_v = obj.deg * N_v;
            dof   = dof_x * dof_v;
            dof_K = obj.deg^2;
            
            dx = ( obj.dim_x.max - obj.dim_x.min ) / N_x;
            dv = ( obj.dim_v.max - obj.dim_v.min ) / N_v;
            
            rho_1_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{1}, [ obj.deg, N_x ] );
            rho_2_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{2}, [ obj.deg, N_x ] );
            rho_3_q = ( obj.Phi / sqrt( dx ) ) * reshape( Q{3}, [ obj.deg, N_x ] );
            
            [ D_q, U_q, T_q ] = Primitive( rho_1_q, rho_2_q, rho_3_q, obj.T_min );
            
            rhs_vDotGradMaxwellian = zeros( dof, 1 );
            
            %% 2D quadrature Weights
            
            w_q = kron( obj.quad_w, obj.quad_w ); % Move to initialization?
            
            %% Defivative of 2D Basis Functions Evaluated in Quadrature Points on Local Element
            
            dPhi_q = kron( obj.dPhi, obj.Phi ) / sqrt( dx * dv ) / ( dx / 2 ); % Move to initialization?
            
            % --- Volume Term ---
            
            i_K = 0;
            for i_x = 1 : N_x
            for i_v = 1 : N_v
                
                i_K = i_K + 1;

                v_q = obj.dim_v.min + ((i_v-1)+.5*(1.+obj.quad_x)) * dv;
                
                M_q = Maxwellian_q( D_q(:,i_x), U_q(:,i_x), T_q(:,i_x), v_q );
                
                v_q = kron( ones(obj.n_quad,1), v_q );
                
                rhs_vDotGradMaxwellian((i_K-1)*dof_K+1:i_K*dof_K)...
                    = 0.25 * dx * dv * ( w_q .* v_q .* M_q )' * dPhi_q;
                
            end
            end
            
            % --- Numerical Fluxes ---
            
            D_L = zeros(N_x+1,1);
            U_L = zeros(N_x+1,1);
            T_L = zeros(N_x+1,1);
            
            D_R = zeros(N_x+1,1);
            U_R = zeros(N_x+1,1);
            T_R = zeros(N_x+1,1);
            
            first = 1:obj.deg;
            last  = (N_x-1)*obj.deg+1:N_x*obj.deg;
            for i_x = 1 : N_x + 1
                
                left  = GetIndices_L( i_x, obj.deg, N_x, obj.BCs_X, first, last );
                right = GetIndices_R( i_x, obj.deg, N_x, obj.BCs_X, first, last );
                
                rho_1_L = ( obj.PhiUp / sqrt( dx ) ) * Q{1}(left);
                rho_2_L = ( obj.PhiUp / sqrt( dx ) ) * Q{2}(left);
                rho_3_L = ( obj.PhiUp / sqrt( dx ) ) * Q{3}(left);
                
                [ D_L(i_x), U_L(i_x), T_L(i_x) ]...
                    = Primitive( rho_1_L, rho_2_L, rho_3_L, obj.T_min );
                
                rho_1_R = ( obj.PhiDn / sqrt( dx ) ) * Q{1}(right);
                rho_2_R = ( obj.PhiDn / sqrt( dx ) ) * Q{2}(right);
                rho_3_R = ( obj.PhiDn / sqrt( dx ) ) * Q{3}(right);
                
                [ D_R(i_x), U_R(i_x), T_R(i_x) ]...
                    = Primitive( rho_1_R, rho_2_R, rho_3_R, obj.T_min );
                
            end
            
            F_num = zeros(obj.n_quad,N_v,N_x+1);
            
            for i_x = 1 : N_x + 1
            for i_v = 1 : N_v
                
                v_q = obj.dim_v.min + ((i_v-1)+.5*(1.+obj.quad_x)) * dv;
                
                if    ( all( v_q <  0 ) )
                    
                    M_q = Maxwellian_q( D_R(i_x), U_R(i_x), T_R(i_x), v_q );
                    
                elseif( all( v_q >= 0 ) )
                    
                    M_q = Maxwellian_q( D_L(i_x), U_L(i_x), T_L(i_x), v_q );
                    
                else
                    assert( false, 'all v_q must be of same sign' )
                end
                
                F_num(:,i_v,i_x) = v_q .* M_q;
                
            end
            end
            
            % --- Surface Term ---
            
            PhiDn_q = kron( obj.PhiDn, obj.Phi ) / sqrt( dx * dv );
            PhiUp_q = kron( obj.PhiUp, obj.Phi ) / sqrt( dx * dv );
            
            i_K = 0;
            for i_x = 1 : N_x
            for i_v = 1 : N_v
                
                i_K = i_K + 1;
                
                rhs_vDotGradMaxwellian((i_K-1)*dof_K+1:i_K*dof_K)...
                    = rhs_vDotGradMaxwellian((i_K-1)*dof_K+1:i_K*dof_K)...
                    - 0.5 * dv * (   PhiUp_q' * ( obj.quad_w .* F_num(:,i_v,i_x+1) )...
                                   - PhiDn_q' * ( obj.quad_w .* F_num(:,i_v,i_x  ) ) );
                
            end
            end
            
            toc
            
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

function [ D, U, T ] = Primitive( rho_1, rho_2, rho_3, T_min )

  D = rho_1;
  U = rho_2 ./ rho_1;
  T = max( ( 2.0 * rho_3 - rho_2.^2 ./ rho_1 ) ./ ( rho_1 ), T_min );

end

function [ M_q ] = Maxwellian_q( D_qx, U_qx, T_qx, v_qv )

    D_q = kron( D_qx, ones(size(v_qv)) );
    U_q = kron( U_qx, ones(size(v_qv)) );
    T_q = kron( T_qx, ones(size(v_qv)) );
    v_q = kron( ones(size(D_qx)), v_qv );
    
    M_q = ( D_q ./ sqrt( 2.*pi.*T_q ) ) .* exp( - (v_q-U_q).^2 ./ ( 2.*T_q ) );

end
