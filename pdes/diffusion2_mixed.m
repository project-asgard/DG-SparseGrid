function pde = diffusion2_mixed(opts)
% Example PDE using the 2D diffusion with anisotropic diffusion tensor. 
% This example PDE is
% time dependent (although not all the terms are time dependent). This
% implies the need for an initial condition. 
% PDE:
% 
% df/dt = div(B*B'\grad f)
%
% where
%
% D = [2 1;1 2];
%
% D = B*B' where B = sqrt(2)/2*[-1 sqrt(3);1 sqrt(3)];
%
% Domain is [0,1]x[0,1]
% Dirichlet boundary condition 
%
% Diffusion terms are dealt with via LDG, i.e., splitting into two first
% order equations:
%
% div(B*B'\grad f)
%
% df/dt = div(B*q) with free boundary BC
%
% and
%
% q=B'\grad f with Dirichlet BCs specified by analytic solution.
%
% Note the choice of mixed variable will perserve a symmetric positive semi
% definite stiffness matrix.  
%
% Run with
%
% implicit
% asgard(@diffusion2_mixed,'timestep_method','BE','dt',0.001,'num_steps',20)

%% Define the dimensions

dim_y = DIMENSION(0,1);
dim_x = DIMENSION(0,1);
dim_x.moment_dV = @(x,p,t,d) 0.*x+1;
dim_y.moment_dV = @(x,p,t,d) 0.*x+1;
dimensions = {dim_x,dim_y};
num_dims = numel(dimensions);

B = sqrt(2)/2*[-1 sqrt(3);1 sqrt(3)];
%B = eye(2);

%% Define the analytic solution (optional).

soln_x = @(x,p,t) sin(pi*x);
soln_y = @(y,p,t) sin(pi*y);
soln_t = @(t,p) cos(t);
soln1 = new_md_func(num_dims,{soln_x,soln_y,soln_t});

solutions = {soln1};

%% Define the boundary conditions

BCL = soln1;
BCR = soln1;

%% Initial conditions

initial_conditions = {soln1};

%% Define the terms of the PDE
%
% Here we have 1 term, with each term having nDims (x and y) operators.

%Surface jacobian matches volume jacobian
dV = @(x,p,t,dat) 0*x+1;

%Set local LF values for 
LF_DIV = 0;
LF_GRD = 0;

%% term1

g1 = @(x,p,t,dat) x.*0+B(1,1);
g2 = @(x,p,t,dat) x.*0+B(1,1);

pterm1 =  DIV(num_dims,g1,'',LF_DIV,'N','N','','','',dV);
pterm2 = GRAD(num_dims,g2,'',LF_GRD,'D','D',BCL,BCR,'',dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = MASS(g1,[],[],dV);
pterm2 = MASS(g2,[],[],dV);
term_y = SD_TERM({pterm1,pterm2});

term1  = MD_TERM(num_dims,{term_x,term_y});

%% term2

g1 = @(x,p,t,dat) x.*0+B(1,1);
g2 = @(x,p,t,dat) x.*0+1;

pterm1 =  DIV(num_dims,g1,'',LF_DIV,'N','N','','','',dV);
pterm2 = MASS(g2,[],[],dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+B(1,2);

pterm1 = MASS(g1,[],[],dV);
pterm2 = GRAD(num_dims,g2,'',LF_GRD,'D','D',BCL,BCR,'',dV);
term_y = SD_TERM({pterm1,pterm2});

term2  = MD_TERM(num_dims,{term_x,term_y});

%% term3

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+B(1,1);

pterm1 = MASS(g1,[],[],dV);
pterm2 = GRAD(num_dims,g2,'',LF_DIV,'D','D',BCL,BCR,'',dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+B(1,2);
g2 = @(x,p,t,dat) x.*0+1;

pterm1 =  DIV(num_dims,g1,'',LF_GRD,'N','N','','','',dV);
pterm2 = MASS(g2,[],[],dV);
term_y = SD_TERM({pterm1,pterm2});

term3  = MD_TERM(num_dims,{term_x,term_y});

%% term4

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = MASS(g1,[],[],dV);
pterm2 = MASS(g2,[],[],dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+B(1,2);
g2 = @(x,p,t,dat) x.*0+B(1,2);

pterm1 =  DIV(num_dims,g1,'',LF_DIV,'N','N','','','',dV);
pterm2 = GRAD(num_dims,g2,'',LF_GRD,'D','D',BCL,BCR,'',dV);
term_y = SD_TERM({pterm1,pterm2});

term4  = MD_TERM(num_dims,{term_x,term_y});

%% term5

g1 = @(x,p,t,dat) x.*0+B(2,1);
g2 = @(x,p,t,dat) x.*0+B(2,1);

pterm1 =  DIV(num_dims,g1,'',LF_DIV,'N','N','','','',dV);
pterm2 = GRAD(num_dims,g2,'',LF_GRD,'D','D',BCL,BCR,'',dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = MASS(g1,[],[],dV);
pterm2 = MASS(g2,[],[],dV);
term_y = SD_TERM({pterm1,pterm2});

term5  = MD_TERM(num_dims,{term_x,term_y});

%% term6

g1 = @(x,p,t,dat) x.*0+B(2,1);
g2 = @(x,p,t,dat) x.*0+1;

pterm1 =  DIV(num_dims,g1,'',LF_DIV,'N','N','','','',dV);
pterm2 = MASS(g2,[],[],dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+B(2,2);

pterm1 = MASS(g1,[],[],dV);
pterm2 = GRAD(num_dims,g2,'',LF_GRD,'D','D',BCL,BCR,'',dV);
term_y = SD_TERM({pterm1,pterm2});

term6  = MD_TERM(num_dims,{term_x,term_y});

%% term7

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+B(2,1);

pterm1 = MASS(g1,[],[],dV);
pterm2 = GRAD(num_dims,g2,'',LF_GRD,'D','D',BCL,BCR,'',dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+B(2,2);
g2 = @(x,p,t,dat) x.*0+1;

pterm1 =  DIV(num_dims,g1,'',LF_DIV,'N','N','','','',dV);
pterm2 = MASS(g2,[],[],dV);
term_y = SD_TERM({pterm1,pterm2});

term7  = MD_TERM(num_dims,{term_x,term_y});

%% term8

g1 = @(x,p,t,dat) x.*0+1;
g2 = @(x,p,t,dat) x.*0+1;

pterm1 = MASS(g1,[],[],dV);
pterm2 = MASS(g2,[],[],dV);
term_x = SD_TERM({pterm1,pterm2});

g1 = @(x,p,t,dat) x.*0+B(2,2);
g2 = @(x,p,t,dat) x.*0+B(2,2);

pterm1 =  DIV(num_dims,g1,'',LF_DIV,'N','N','','','',dV);
pterm2 = GRAD(num_dims,g2,'',LF_GRD,'D','D',BCL,BCR,'',dV);
term_y = SD_TERM({pterm1,pterm2});

term8  = MD_TERM(num_dims,{term_x,term_y});

%% Penalty terms

%% term9

g1 = @(x,p,t,dat) x.*0+1;
pterm1 =  DIV(num_dims,g1,'',-1,'N','N','','','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*0+1;
pterm1 = MASS(g1,[],[],dV);
term_y = SD_TERM({pterm1});

term9  = MD_TERM(num_dims,{term_x,term_y});

%% term10

g1 = @(x,p,t,dat) x.*0-1;
pterm1 =  DIV(num_dims,g1,'',0,'D','D','','','',dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*0+1;
pterm1 = MASS(g1,[],[],dV);
term_y = SD_TERM({pterm1});

term10  = MD_TERM(num_dims,{term_x,term_y});

%% term11

g1 = @(x,p,t,dat) x.*0+1;
pterm1 = MASS(g1,[],[],dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*0+1;
pterm1 =  DIV(num_dims,g1,'',-1,'N','N','','','',dV);
term_y = SD_TERM({pterm1});

term11  = MD_TERM(num_dims,{term_x,term_y});

%% term12

g1 = @(x,p,t,dat) x.*0+1;
pterm1 = MASS(g1,[],[],dV);
term_x = SD_TERM({pterm1});

g1 = @(x,p,t,dat) x.*0-1;
pterm1 =  DIV(num_dims,g1,'',0,'D','D','','','',dV);
term_y = SD_TERM({pterm1});

term12  = MD_TERM(num_dims,{term_x,term_y});

%%%Without penalty
terms = {term1, term2, term3, term4, term5, term6, term7, term8};

%%%With penalty
%terms = {term1, term2, term3, term4, term5, term6, term7, term8, ...
%         term9, term10, term11, term12};

%% Define some parameters and add to pde object.
%  These might be used within the various functions below.

params.parameter1 = 0;
params.parameter2 = 1;

%% Define sources

sc_x = @(x,p,t,d) sin(pi*x);
sc_y = @(x,p,t,d) sin(pi*x);
sc_t = @(t,p) -sin(t);
s1 = new_md_func(num_dims,{sc_x,sc_y,sc_t});

sc_x = @(x,p,t,d) sin(pi*x);
sc_y = @(x,p,t,d) sin(pi*x);
sc_t = @(t,p) 4*pi^2*cos(t);
s2 = new_md_func(num_dims,{sc_x,sc_y,sc_t});

sc_x = @(x,p,t,d) cos(pi*x);
sc_y = @(x,p,t,d) cos(pi*x);
sc_t = @(t,p) -2*pi^2*cos(t);
s3 = new_md_func(num_dims,{sc_x,sc_y,sc_t});

sources = {s1,s2,s3};

%% Define a function to set dt

    function dt=set_dt(pde,CFL)
        dims = pde.dimensions;      
        % for Diffusion equation: dt = C * dx^2
        lev = dims{1}.lev;
        dx = 1/2^lev;
        dt = CFL*dx^2;
    end

%% Construct PDE

pde = PDE(opts,dimensions,terms,[],sources,params,@set_dt,[],initial_conditions,solutions);

end

