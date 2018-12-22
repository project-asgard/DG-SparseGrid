function pde=Vlasov8
% Numerical Example for Vlasov Equation
% This test has given E and non-zero source term

% parameters

pde.dims = ['x','v'];
pde.BCL = [0,0];
pde.BCR = [0,0];
pde.domainMin = [-1,-10];
pde.domainMax = [+1,+10];
pde.lev = [2,2];

%% 
% Setup the v.d_dx (v.MassV . GradX)
% We construct in the order specified by p.dims, i.e., 
% [GradX, v.MassV]
term2.name = 'v.d_dx';
term2.type = [1,2]; % [grad,mass]
term2.g = {@(x,t,dat) 1, @(v,t,dat) v};
term2.TD = [0,0];
term2.dat = {[],[]} % These are to be filled within the workflow for now
term2.LF = 0;

%% 
% Setup the E.d_dv (E.MassX . GradV)
% We construct in the order specified by p.dims.
term3.name = 'E.d_dv';
term3.type = [2,1]; % [grad,mass]
term3.g = {@(x,t,dat) dat, @(v,t,dat) 1};
term3.TD = [1,0];
term3.dat = {[],[]}; % These are to be filled within the workflow for now
term3.LF = 0;


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
