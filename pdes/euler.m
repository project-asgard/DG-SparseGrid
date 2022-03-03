function pde = euler(opts)

%% Define the dimensions

dim_x = DIMENSION(0.0,1.0);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

%% Initial Conditions
% case = 1 % Advection of sinusoidal density profile
% case = 2 % Linear sound waves
% case = 3 % Sod's Riemann problem

switch opts.case_
    case 1
        
        U_0 = 1.0e-0;
        T_0 = 1.0e-10;

        % rho_1:
        ic_x = @(x,p,t,d) 0.5*sin(2*pi.*x)+1.0;
        ic_t = @(t,p) 0*t+1; 
        ic_rho_1 = {new_md_func(num_dims,{ic_x,ic_t})};

        % rho_2:
        ic_x = @(x,p,t,d) (0.5*sin(2*pi.*x)+1.0).*U_0;
        ic_t = @(t,p) 0*t+1; 
        ic_rho_2 = {new_md_func(num_dims,{ic_x,ic_t})};

        % rho_3:
        ic_x = @(x,p,t,d) 0.5*(0.5*sin(2*pi.*x)+1.0).*(U_0^2+T_0);
        ic_t = @(t,p) 0*t+1; 
        ic_rho_3 = {new_md_func(num_dims,{ic_x,ic_t})};
        
    case 2
        
        N_0 = 1.0;
        T_0 = 1.0/3.0;
        C_0 = sqrt(3.0*T_0);
        Amp = 1.0e-6;
        
        % rho_1:
        ic_x = @(x,p,t,d) N_0+Amp*sin(2*pi.*x)/C_0^2;
        ic_t = @(t,p) 0*t+1; 
        ic_rho_1 = {new_md_func(num_dims,{ic_x,ic_t})};

        % rho_2:
        ic_x = @(x,p,t,d) Amp*sin(2*pi.*x)/C_0;
        ic_t = @(t,p) 0*t+1; 
        ic_rho_2 = {new_md_func(num_dims,{ic_x,ic_t})};

        % rho_3:
        ic_x = @(x,p,t,d) 0.5*(N_0*T_0+Amp*sin(2*pi.*x));
        ic_t = @(t,p) 0*t+1; 
        ic_rho_3 = {new_md_func(num_dims,{ic_x,ic_t})};
        
    case 3
        
        x_D = 64/128;
        
        N_L = 1.0;
        U_L = 0.0;
        T_L = 1.0;
        
        rho_1_L = N_L;
        rho_2_L = N_L * U_L;
        rho_3_L = 0.5 * N_L * (U_L^2+T_L);
        
        N_R = 0.125;
        U_R = 0.0;
        T_R = 0.8;
        
        rho_1_R = N_R;
        rho_2_R = N_R * U_R;
        rho_3_R = 0.5 * N_R * (U_R^2+T_R);
        
        % rho_1:
        ic_x = @(x,p,t,d) rho_1_R + (rho_1_L-rho_1_R).*(x <= x_D);
        ic_t = @(t,p) 0*t+1; 
        ic_rho_1 = {new_md_func(num_dims,{ic_x,ic_t})};
        
        % rho_2:
        ic_x = @(x,p,t,d) rho_2_R + (rho_2_L-rho_2_R).*(x <= x_D);
        ic_t = @(t,p) 0*t+1; 
        ic_rho_2 = {new_md_func(num_dims,{ic_x,ic_t})};

        % rho_3:
        ic_x = @(x,p,t,d) rho_3_R + (rho_3_L-rho_3_R).*(x <= x_D);
        ic_t = @(t,p) 0*t+1; 
        ic_rho_3 = {new_md_func(num_dims,{ic_x,ic_t})};
        
    otherwise
        
end

initial_conditions = {ic_rho_1,ic_rho_2,ic_rho_3};

%% Construct moments

%mass moment
moment_func = new_md_func(num_dims,{@(x,p,t) 0*x+1,@(p,t)   0*t+1});
moment0 = MOMENT({moment_func});

moments = {moment0};

lev = opts.lev;
deg = opts.deg;
h = (dim_x.max-dim_x.min)/2^lev;
Jacobi = h/2.0;

[quad_x,quad_w] = lgwt(default_quad_number(deg),-1,1);
p_val = lin_legendre(quad_x,deg) * 1/sqrt(h);

rho_1_initial = @(x) 0*x+1.0;
rho_2_initial = @(x) 0*x+0.0;
rho_3_initial = @(x) 0*x+0.5;

rho_1 = zeros(deg*2^lev,1);
rho_2 = zeros(deg*2^lev,1);
rho_3 = zeros(deg*2^lev,1);

for i = 1 : deg*2^lev
    
    quad_xi = dim_x.min+0.5*(1.0+quad_x)*h;
    
    rho_1((i-1)*deg+1:i*deg) = ( p_val' * ( quad_w .* rho_1_initial(quad_xi) ) ) .* Jacobi;
    rho_2((i-1)*deg+1:i*deg) = ( p_val' * ( quad_w .* rho_2_initial(quad_xi) ) ) .* Jacobi;
    rho_3((i-1)*deg+1:i*deg) = ( p_val' * ( quad_w .* rho_3_initial(quad_xi) ) ) .* Jacobi;
    
end
rho = {rho_1,rho_2,rho_3};

nlinterms = {@(pde,opts,rho) euler_dg_eval_rhs(pde,opts,rho)};

%% Define a function to set dt

function dt=set_dt(pde,CFL)
    dims = pde.dimensions;      
    % for hyperbolic equation: dt = C * dx
    lev = dims{1}.lev;
    dx = 1/2^lev;
    dt = CFL*dx;
end

pde = PDE(opts,dimensions,[],nlinterms,[],[],[],@set_dt,[],initial_conditions,[],moments);

end

