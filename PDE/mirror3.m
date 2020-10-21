function pde = mirror3(opts)
% Three-dimensional magnetic mirror from the FP paper - evolution of the ion velocity and spatial dependence
% of f in the presence of Coulomb collisions with background electrons
% 
% df/dt == -dB/ds f - vcos(z)df/dz + nu_D/(2*sin(z)) d/dz ( sin(z) df/dz ) + 1/v^2 (d/dv(flux_v))
%a
% flux_v == v^3[(m_a/(m_a + m_b))nu_s f) + 0.5*nu_par*v*d/dv(f)]
%
% Run with
%
% asgard(@mirror3,'timestep_method','BE','case',3)

% choose a case
%
% case 1 - maxwellian offset and different Temp.
% case 2 - max. no offset and different Temp.
% case 3 - max. offset and same Temp.

params = mirror_parameters();

switch opts.case_
    case 1
        params.a.T_eV = 0.05*params.b.T_eV; %Target temperature in Kelvin
        params.init_cond_v = @(v) params.maxwell(v,params.v_th(params.a.T_eV,params.a.m),1e6);
    case 2
        params.a.T_eV = 0.05*params.b.T_eV;
        params.init_cond_v = @(v) params.maxwell(v,0,1e6);
    case 3
        params.a.T_eV = 1e3;
        params.init_cond_v = @(v) v.*0 + 1;%params.maxwell(v,params.v_th(params.a.T_eV,params.a.m),1e6);
end

%v_ = 10.^[-1:.1:7];
%loglog(0.5*a.m*v_.^2/e,nu_D(v_,a,b))
%hold on
%loglog(0.5*a.m*v_.^2/e,nu_s(v_,a,b))
%loglog(0.5*a.m*v_.^2/e,nu_par(v_,a,b))
%xlim([0.1,1e2]);
%ylim([1e4,1e11]);

% Domain is (a,b)

% The function is defined for the plane
% x = a and x = b
BCL_fList = { ...
    @(v,p,t) v.*0, ... 
    @(z,p,t) p.init_cond_z(z), ...
    @(s,p,t) p.init_cond_s(s), ...
    @(t,p)   t.*0 + 1
    };

BCR_fList = { ...
    @(v,p,t) p.init_cond_v(v), ... % 
    @(z,p,t) p.init_cond_z(z), ...
    @(s,p,t) p.init_cond_s(s), ...
    @(t,p) t.*0 + 1
    };

%E = 1.0; %parallel Electric field

%% Define the dimensions

dim_v.name = 'v';
dim_v = DIMENSION(0,3e7);
dim_v.init_cond_fn = @(v,p,t) p.init_cond_v(v);
dim_v.jacobian = @(v,p,t) 2.*pi.*v.^2;

dim_z = DIMENSION(-pi,+pi);
dim_z.name = 'z';
dim_z.init_cond_fn = @(z,p,t) p.init_cond_z(z);
dim_z.jacobian = @(z,p,t) sin(z);

dim_s = DIMENSION(0,5);
dim_s.init_cond_fn = @(s,p,t) p.init_cond_s(s);
dim_s.jacobian = @(s,p,t) s.*0 + 1;

dimensions = {dim_v,dim_z,dim_s};
num_dims = numel(dimensions);

%% Define the terms of the PDE

%% Advection  term
% termS1 == -vcos(z)*df/ds
% termS1 == q(v)*r(z)*w(s)
% q(v) == g1(v)  [mass, g1(p) = v,  BC N/A]
% r(z) == g2(z) [mass, g2(z) = cos(z),  BC N/A]
% w(s) == d/ds g3(s) f [grad, g3(s) = -1, BCL= D, BCR=D]
g1 = @(v,p,t,dat) v;
g2 = @(z,p,t,dat) cos(z);
g3 = @(s,p,t,dat) s.*0 - 1;

pterm1 = MASS(g1);
pterm2 = MASS(g2);
pterm3 = GRAD(num_dims,g3,-1,'D','D', BCL_fList, BCR_fList);

term1_s = TERM_1D({pterm1,pterm2,pterm3});
termS1   = TERM_ND(num_dims,{term1_s,[],[]});

%% Mass term
% termS2 == -vcos(z)dB/ds f
% termS1 == q(v)*r(z)*w(s)
% q(v) == g1(v)  [mass, g1(p) = v,  BC N/A]
% r(z) == g2(z) [mass, g2(z) = cos(z),  BC N/A]
% w(s) == g3(s) f [mass, g3(s) = -dB/ds/B, BCL= D, BCR=D]

g1 = @(v,p,t,dat) v;
g2 = @(z,p,t,dat) cos(z);
g3 = @(s,p,t,dat) -p.dB_ds(s)./p.B_func(s);

pterm1 = MASS(g1);
pterm2 = MASS(g2);
pterm3 = MASS(g3);

termB1 = TERM_1D({pterm1,pterm2,pterm3});

termS2 = TERM_ND(num_dims,{termB1,[],[]});

%% 
% termC == nu_D/(2*sin(z))*d/dz sin(z)*df/dz
%
% becomes 
%
% termC == g1(v) g2(z) q(z)   [mass, g1(p) = nu_D(v), g2(z) = 1/(2sin(z))  BC N/A]
%   q(z) == d/dz g3(z) r(z)   [grad, g3(z) =  sin(z), BCL=D,BCR=D]
%   r(z) == d/dp g4(z) f(z)   [grad, g3(p) = 1,      BCL=N,BCR=N]

g1 = @(v,p,t,dat) p.nu_D(v,p.a,p.b);
g2 = @(z,p,t,dat) 1./(2*sin(z));
g3 = @(z,p,t,dat) sin(z);
g4 = @(z,p,t,dat) z.*0 + 1;
pterm1  = MASS(g1);
pterm2  = MASS(g2);
pterm3  = GRAD(num_dims,g3,+1,'D','D');
pterm4  = GRAD(num_dims,g4,-1,'N', 'N');
termC_z = TERM_1D({pterm1,pterm2,pterm3,pterm4});
termC   = TERM_ND(num_dims,{termC_z,[],[]});

% term V1 == 1/v^2 d/dv(v^3(m_a/(m_a + m_b))nu_s f))
% term V1 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)f(v))   [grad, g2(v) = v^3(m_a/(m_a + m_b))nu_s, BCL= N, BCR=D]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^3*p.a.m.*p.nu_s(v,p.a,p.b)./(p.a.m + p.b.m);

pterm1  = MASS(g1);
pterm2  = GRAD(num_dims,g2,-1,'N','D', BCL_fList, BCR_fList);
termV_s = TERM_1D({pterm1,pterm2});
termV1  = TERM_ND(num_dims,{termV_s,[],[]});

%%
% term V2 == 1/v^2 d/dv(v^4*0.5*nu_par*d/dv(f))
% term V2 == g(v) q(v)      [mass, g(v) = 1/v^2,  BC N/A]
% q(v) == d/dv(g2(v)r(v))   [grad, g2(v) = v^4*0.5*nu_par, BCL= D, BCR=D]
% r(v) = d/dv(g3(v)f)       [grad, g3(v) = 1, BCL=N, BCR=N]

g1 = @(v,p,t,dat) 1./v.^2;
g2 = @(v,p,t,dat) v.^4.*0.5.*p.nu_par(v,p.a,p.b);
g3 = @(v,p,t,dat) v.*0 + 1;

pterm1      = MASS(g1);
pterm2      = GRAD(num_dims,g2,-1,'D','N');
pterm3      = GRAD(num_dims,g3,+1,'N','D', BCL_fList, BCR_fList);
termV_par   = TERM_1D({pterm1,pterm2,pterm3});
termV2      = TERM_ND(num_dims,{termV_par,[],[]});

terms = {termV1,termV2,termC,termS1};

%% Define sources

sources = {};

%% Define the analytic solution (optional).
% This requires nDims+time function handles.

analytic_solutions_1D = { ...    
    @(v,p,t) p.init_cond_v(v), ...
    @(z,p,t) p.init_cond_z(z), ...
    @(s,p,t) p.analytic_solution_s(s,p,t), ...
    @(t,p) t.*0 + 1; %pitch_t(t)
    };

%% Define function to set time step
    function dt=set_dt(pde,CFL)      
        Lmax = pde.dimensions{1}.max;
        LevX = pde.dimensions{1}.lev;
        dt = Lmax/2^LevX*CFL;
    end


%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,analytic_solutions_1D);

end
