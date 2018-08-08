function pde=Vlasov7
% Numerical Example for Vlasov Equation
%--------------------------------------------------
% d/dt(f) + v*d/dx(f) + E*d/dv(f)=0
% term1 denotes v*d/dx(f)
%   term1 = [formT1,formT2]
%   formT1 = FuncMass :: \int_v v*phi_i*phi_j dv
%   formT2 = FuncGrad :: \int_x phi_i'*phi_j dx
% term2 denotes E*d/dv(f)
%   term2 = [formT1,formT2]
%   formT1 = FuncGrad :: \int_v phi_i'*phi_j dv
%   formT2 = FuncMass :: \int_x E*phi_i*phi_j dx
% Note: All terms will be constructed by weak forms
%--------------------------------------------------
Dim = 2;

formT1.dim = 1;
formT1.type = 'FuncMass';
formT1.G = @(x)x;

formT2.dim = 1;
formT2.type = 'FuncGrad';
formT2.G = @(x)1;

term1 = {formT1,formT2};

formT1.dim = 1;
formT1.type = 'FuncMass';
formT1.G = @(x)2*x;

formT2.dim = 1;
formT2.type = 'FuncGrad';
formT2.G = @(x)2;

term2 = {formT1,formT2};


terms = {term1,term2};

pde.terms = terms;


% This test has given E and non-zero source term


% parameters
k_0=0.5;
A=1;
Vmax=5;Lmax=1;

params.k_0 = k_0;
params.A = A;
params.Lmax = Lmax;
params.Vmax = Vmax;

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

pde.checkAnalytic = 1;

pde.exactE = @exactE;

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

end

function f=Fx_0(x, params)
% Initial condition for x variable

A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=x.*(Lmax-x);
end

function f=Fv_0(v, params)
% Initial condition for v variable

A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=v-v;%+1;%(v+Vmax).*(v-Vmax);
end
function f=Fxv_0(x,v, params)
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=Fv_0(v).*Fx_0(x);
end
function f=Ex(x,  params)
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=x-x+1;
end
function f=Et(t, params)
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=t-t+1;
end
function f=E(x,t, params)
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=Ex(x).*Et(t);
end
function f=F(x,v,t,  params)
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=t.*x.*(1-x).*(v-Vmax).*(v+Vmax);%.*x.*(Lmax-x);%t*(v-Vmax).*(v+Vmax).*x.*(Lmax-x);
end
function f=rho(x,t, params)
% we do note use rho in this test, so I just write an arbitrary function
A = params.A;
k_0 = params.k_0;
Lmax = params.Lmax;
Vmax = params.Vmax;

f=x-x+1;
end

function f = exactE(x)
% Exact solution for E
f=x-x+1;
end
% source term--fully seperable functions
% source = source1+source2+source3
% source term 1
function f = source1x(x)
Vmax=5;Lmax=1;
f = x.*(1-x);
end
function f = source1t(t)
Vmax=5;Lmax=1;
f = t-t+1;%+1;
end
function f = source1v(v)
Vmax=5;Lmax=1;
f = (v-Vmax).*(v+Vmax);
end
function f = source1(x,v,t)
Vmax=5;Lmax=1;
f = source1x(x).*source1v(v).*source1t(t);
end

% source term 2
function f = source2x(x)
Vmax=5;Lmax=1;
f = (1-2*x);
end
function f = source2t(t)
Vmax=5;Lmax=1;
f = t;%-t+1;
end
function f = source2v(v)
Vmax=5;Lmax=1;
f = v.*(v-Vmax).*(v+Vmax);
end
function f = source2(x,v,t)
Vmax=5;Lmax=1;
f = source2x(x).*source2v(v).*source2t(t);
end

% source term 3
function f = source3x(x)
Vmax=5;Lmax=1;
f = x.*(1-x);%.*exactE(x);
end
function f = source3t(t)
Vmax=5;Lmax=1;
f = t;
end
function f = source3v(v)
Vmax=5;Lmax=1;
f = 2*v;
end
function f = source3(x,v,t)
Vmax=5;Lmax=1;
f = source3x(x).*source3v(v).*source3t(t);
end

% Exact F
function f=ExactFx(x)
Vmax=5;Lmax=1;
f = x.*(1-x);
end
function f=ExactFv(v)
Vmax=5;Lmax=1;
f = (v-Vmax).*(v+Vmax);
end
function f=ExactFt(t)
Vmax=5;Lmax=1;
f=t;%+1;
end
function f=ExactF(x,v,t)
Vmax=5;Lmax=1;
f = ExactFx(x).*ExactFv(v).*ExactFt(t);
end
