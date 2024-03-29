function pde = fokkerplanck1_momentum_R(opts)
% Combining momentum and pitch angle dynamics for the E term
%
% d/dt f(p,z) == -div(flux_R)
%
% flux_R is the flux due to radiative damping
%
% -div(flux_R) == termR1
%
% setting E=2 and z=0.5, and expanding the d/dz term via product rule, then setting d/dz=0,
% we get for the momentum direction ...
%
% termR1 == -1/p^2 (d/dp p^2 3/4*f(p)) [pterm1: div (g(p)=-3/4,-1,BCL=N,BCR=N,dV=p^2]
%
% p is a spherical coordinate (p,th,ph), so the div and grad look like
%
% div[] = 1/p^2 * d/dp * p^2[], grad[] = d/dp[]
%
% and the volument_element dV = p^2
%
%
% Run with
%
% asgard(@fokkerplanck1_momentum_R,'timestep_method','matrix_exponential','dt',0.1,'num_steps',10,'case',1,'lev',4,'deg',4)
%
% case = 1 % offset maxwellian
%
% Notes:
% DLG - at present this is a non-conservative PDE, so the mass will
% increase. I may fix this in future as I did with the
% fokkerplanck1_momentum_E_div test. Updated now. 

params = fokkerplanck_parameters(opts);

%% Setup the dimensions

dV_p = @(x,p,t,d) x.^2;
dim_p = DIMENSION(0,+10);
dim_p.moment_dV = dV_p;
dimensions = {dim_p};
num_dims = numel(dimensions);

%% Define the analytic solution (optional)

    function ans = soln(x,t)
        A = exp(-(-4+1./4*(-3.*t+4.*x)).^2);
        y = -3.*t+4.*x;
%         B = sign(y).*abs(y).^(5./3);
%         C = (8.*2.^(1./3).*x.^(5./3));
        D = 16.*x.^2;       
%         ans = A .* B ./ C;
        ans = A .* y ./ D;
    end

switch opts.case_
    case 1
        soln_p = @(x,p,t) soln(x,t);
end

soln1 = new_md_func(num_dims,{soln_p});
solutions = {soln1};

%% Define the initial conditions

switch opts.case_
    case 1
        initial_conditions = {soln1};
end

%% Define the terms of the PDE

%% -div(flux_R) == termR1

% termR1 == -1/p^2 (d/dp p^2 3/4*f(p)) [div, g1(p) = -3/4,   BCL=N,BCR=N]

g1     = @(x,p,t,d) x.*0-3/4;
pterm1 = DIV (num_dims,g1,'',-1,'N','N',soln1,soln1,'',dV_p);

termR1_p = SD_TERM({pterm1});
termR1   = MD_TERM(num_dims,{termR1_p});

terms = {termR1};

%% Define sources

sources = {};

%% Define function to set time step
    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;
        xRange = dims{1}.max-dims{1}.min;
        lev = dims{1}.lev;
        dx = xRange/2^lev;
        dt = CFL * dx;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end


