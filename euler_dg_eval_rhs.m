function [ rhs_W ] = euler_dg_eval_rhs( pde, opts, rho_W )

dim = pde.dimensions{1};

d_v = 1.0; % Velocity space dimensions hard coded.

lev = opts.lev;
deg = opts.deg;

N = 2^lev;
h = (dim.max-dim.min) / N;
dof_1D = deg * N;
Jacobi = h/2;

[quad_x,quad_w] = lgwt(default_quad_number(deg),-1,1);

%% Basis functions and derivatives evaluated interior to element
%  p_val (:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg

p_val  = lin_legendre (quad_x,deg) * 1/sqrt(h);
dp_val = lin_dlegendre(quad_x,deg) * 1/sqrt(h) / Jacobi;

%% Basis functions eveluated on element edges
%  p_L(:) is 1 by deg
%  p_R(:) is 1 by deg

p_L = lin_legendre(-1,deg) * 1/sqrt(h);
p_R = lin_legendre(+1,deg) * 1/sqrt(h);

%% Convert from wavelet space to real space

FMWT = OperatorTwoScale_wavelet2(deg,lev);

rho_real = cell(numel(rho_W),1);
for i = 1 : numel(rho_W)
  rho_real{i} = FMWT' * rho_W{i};
end

rhs_1_real = zeros(dof_1D,1);
rhs_2_real = zeros(dof_1D,1);
rhs_3_real = zeros(dof_1D,1);

%% Loop over elements to construct volume term of rhs in real space
for i = 1 : N
    
    %% Index Ranges
    
    current = (i-1)*deg+1:i*deg;
        
    %% Volume Term
    
    rho_1 = p_val * rho_real{1}(current);
    rho_2 = p_val * rho_real{2}(current);
    rho_3 = p_val * rho_real{3}(current);
    
    [ D, U, T ] = Primitive( rho_1, rho_2, rho_3, d_v );
    
    [ F_1, F_2, F_3 ] = FluxVector( D, U, T, d_v );
    
    rhs_1_real(current,1) = ( dp_val' * ( quad_w .* F_1 ) ) .* Jacobi;
    rhs_2_real(current,1) = ( dp_val' * ( quad_w .* F_2 ) ) .* Jacobi;
    rhs_3_real(current,1) = ( dp_val' * ( quad_w .* F_3 ) ) .* Jacobi;
    
end

%% Compute Numerical Fluxes

F_1_Num = zeros(N+1,1);
F_2_Num = zeros(N+1,1);
F_3_Num = zeros(N+1,1);

for i = 1 : N + 1
    
    %% Index Ranges
    
    first = 1:deg;
    last  = (N-1)*deg+1:N*deg;
    
    if i == 1
        left = last;
    else
        left = (i-2)*deg+1:(i-1)*deg;
    end

    if i == N + 1
        right = first;
    else
        right = (i-1)*deg+1:i*deg;
    end
    
    %% Left State
    
    rho_1_L = p_R * rho_real{1}(left);
    rho_2_L = p_R * rho_real{2}(left);
    rho_3_L = p_R * rho_real{3}(left);
    
    [ D_L, U_L, T_L ] = Primitive( rho_1_L, rho_2_L, rho_3_L, d_v );
    
    Cs_L = sqrt( (2.0+d_v) * T_L / d_v );
    
    Lambda_L = max( [ abs(U_L+Cs_L), abs(U_L-Cs_L) ] );
    
    [ F_1_L, F_2_L, F_3_L ] = FluxVector( D_L, U_L, T_L, d_v );
    
    %% Right State
    
    rho_1_R = p_L * rho_real{1}(right);
    rho_2_R = p_L * rho_real{2}(right);
    rho_3_R = p_L * rho_real{3}(right);
    
    [ D_R, U_R, T_R ] = Primitive( rho_1_R, rho_2_R, rho_3_R, d_v );
    
    Cs_R = sqrt( (2.0+d_v) * T_R / d_v );
    
    Lambda_R = max( [ abs(U_R+Cs_R), abs(U_R-Cs_R) ] );
    
    [ F_1_R, F_2_R, F_3_R ] = FluxVector( D_R, U_R, T_R, d_v );
    
    %% Numerical Flux
    
    alpha = max( [ Lambda_L, Lambda_R ] );
    
    F_1_Num(i) = NumericalFlux_LLF( rho_1_L, rho_1_R, F_1_L, F_1_R, alpha );
    F_2_Num(i) = NumericalFlux_LLF( rho_2_L, rho_2_R, F_2_L, F_2_R, alpha );
    F_3_Num(i) = NumericalFlux_LLF( rho_3_L, rho_3_R, F_3_L, F_3_R, alpha );
    
end

%% Loop over elements to construct surface term of rhs in real space
for i = 1 : N
    
    %% Index Ranges
    
    current = (i-1)*deg+1:i*deg;
    
    %% Surface Term
    
    rhs_1_real(current,1)...
      = rhs_1_real(current,1)...
      - ( p_R' .* F_1_Num(i+1) - p_L' .* F_1_Num(i) );
    
    rhs_2_real(current,1)...
      = rhs_2_real(current,1)...
      - ( p_R' .* F_2_Num(i+1) - p_L' .* F_2_Num(i) );
    
    rhs_3_real(current,1)...
      = rhs_3_real(current,1)...
      - ( p_R' .* F_3_Num(i+1) - p_L' .* F_3_Num(i) );
    
end

rhs_real = {rhs_1_real,rhs_2_real,rhs_3_real};


%% Convert from real space to wavelet space

rhs_W = cell(numel(rho_W),1);
for i = 1 : numel(rho_W)
  rhs_W{i} = FMWT * rhs_real{i};
end

end

function [ N, U, T ] = Primitive( rho_1, rho_2, rho_3, d_v )

N = rho_1;
U = rho_2 ./ rho_1;
T = ( 2.0 * rho_3 - rho_2.^2 ./ rho_1 ) ./ ( d_v * rho_1 );

end

function [ F_1, F_2, F_3 ] = FluxVector( N, U, T, d_v )

F_1 = N .* U;
F_2 = N .* ( U.^2 + T );
F_3 = 0.5 .* N .* ( U.^2 + (d_v+2) * T ) .* U;

end

function [ F_Num ] = NumericalFlux_LLF( rho_L, rho_R, F_L, F_R, alpha )

F_Num = 0.5 .* ( F_R + F_L - alpha .* ( rho_R - rho_L ) );

end