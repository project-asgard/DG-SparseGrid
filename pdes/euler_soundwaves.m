function pde = euler_soundwaves(opts)

%% Define the dimensions

dim_x = DIMENSION(0.0,1.0);
dim_x.moment_dV = @(x,p,t,dat) 0*x+1;
dimensions = {dim_x};
num_dims = numel(dimensions);

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

pde = {};
%pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions,moments);

F = euler_dg_eval_F(pde,opts,rho);

end

