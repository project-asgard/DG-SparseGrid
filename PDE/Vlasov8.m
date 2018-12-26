function pde=Vlasov8
% Numerical Example for Vlasov Equation
% This test has given E and non-zero source term

% PDE parameters

dim_x.name = 'x';
dim_x.BCL = 0;
dim_x.BCR = 0;
dim_x.domainMin = -1;
dim_x.domainMax = +1;
dim_x.lev = 2;
dim_x.deg = 2;
dim_x.FMWT = [];

dim_v.name = 'v';
dim_v.BCL = 0;
dim_v.BCR = 0;
dim_v.domainMin = -10;
dim_v.domainMax = +10;
dim_v.lev = 2;
dim_v.deg = 2;
dim_v.FMWT = [];

pde.dimensions = {dim_x, dim_v};

%% 
% Setup the v.d_dx (v.MassV . GradX) term

term2_x.type = 1; % grad
term2_x.G = @(x,t,dat) 1;
term2_x.TD = 0;
term2_x.dat = []; % These are to be filled within the workflow for now
term2_x.LF = 0;

term2_v.type = 2; % mass
term2_v.G = @(v,t,dat) v;
term2_v.TD = 0;
term2_v.dat = []; % These are to be filled within the workflow for now
term2_v.LF = 0;

term2.name = 'v.d_dx';
term2 = {term2_x, term2_v};

%% 
% Setup the E.d_dv (E.MassX . GradV) term

term3_x.type = 2; % mass
term3_x.G = @(x,t,dat) dat;
term3_x.TD = 1;
term3_x.dat = []; % These are to be filled within the workflow for now
term3_x.LF = 0;

term3_v.type = 1; % mass
term3_v.G = @(v,t,dat) 1;
term3_v.TD = 0;
term3_v.dat = []; % These are to be filled within the workflow for now
term3_v.LF = 0;

term3.name = 'E.d_dv';
term3 = {term3_x, term3_v};


%%
% Add terms to PDE

pde.terms = {term2, term3};

p.Lmin=-1;
p.Lmax=+1;
p.Vmin=-10;
p.Vmax=+10;

p.TEND = 4;

params = p;

pde.Fx_0 = @Fx_0;
pde.Fv_0 = @Fv_0;
pde.Fxv_0 = @Fxv_0;
pde.Ex = @Ex;
pde.Et = @Et;
pde.E = @E;
pde.rho = @rho;
pde.params = params;

pde.solvePoisson = 0;
pde.applySpecifiedE = 1;
pde.implicit = 0;
pde.checkAnalytic = 1;

pde.source1x = @source1x;
pde.source1v = @source1v;
pde.source1t = @source1t;
pde.source1 = @source1;

pde.source2x = @source2x;
pde.source2v = @source2v;
pde.source2t = @source2t;
pde.source2 = @source2;

pde.source3x = @source3x;
pde.source3v = @source3v;
pde.source3t = @source3t;
pde.source3 = @source3;

pde.ExactFx = @ExactFx;
pde.ExactFv = @ExactFv;
pde.ExactFt = @ExactFt;
pde.ExactF = @ExactF;

end

function f=Fx_0(x,p)
% Initial condition for x variable
f=x.*0;
end
function f=Fv_0(v,p)
% Initial condition for v variable
f=v.*0;
end
function f=Fxv_0(x,v,p)
f=Fv_0(v).*Fx_0(x);
end

% Apply this specific E field
function f=Ex(x, p)
f=cos(pi*x);
end
function f=Et(t,p)
f=cos(t);
end
function f=E(x,t,p)
f=Ex(x).*Et(t);
end

% source term--fully seperable functions
% source = source1+source2+source3

% source term 1
function f = source1t(t)
f = cos(t);
end
function f = source1x(x,p)
f = sin(pi*x);
end
function f = source1v(v,p)
f = sin(pi*v/5);
end
function f = source1(x,v,t)
f = source1x(x).*source1v(v).*source1t(t);
end

% source term 2
function f = source2t(t)
f = sin(t);
end
function f = source2x(x,p)
f = cos(pi*x);
end
function f = source2v(v,p)
f = pi*v.*sin(pi*v/5);
end
function f = source2(x,v,t)
f = source2x(x).*source2v(v).*source2t(t);
end

% source term 3
function f = source3t(t)
f = cos(t).*sin(t);
end
function f = source3x(x,p)
f = cos(pi*x).*sin(pi*x);
end
function f = source3v(v,p)
f = 1/5*pi*cos(pi*v/5);
end
function f = source3(x,v,t)
f = source3x(x).*source3v(v).*source3t(t);
end


% Exact F
function f=ExactFt(t)
f=sin(t);
end
function f=ExactFx(x,p)
f = sin(pi*x);
end
function f=ExactFv(v,p)
f = sin(pi*v/5);
end
function f=ExactF(x,v,t)
f = ExactFx(x).*ExactFv(v).*ExactFt(t);
end
