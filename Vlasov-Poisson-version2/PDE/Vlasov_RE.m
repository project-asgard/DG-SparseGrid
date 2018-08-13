function pde=Vlasov_RE
% Numerical Example for Vlasov Equation
%--------------------------------------------------
% Original equation: d/dt(f) + v*d/dx(f) + E*d/dv(f)=0
% This function's equation: d/dt(f) + EComp + RadComp = CollComp
% v replaced with p (momentum)
% x replaced with xi (angle)
% FuncMass(G) :: \int_x G*phi_i*phi_j dx
% FuncGrad(G) :: \int_x G*phi_i'*phi_j dp
% FuncGrad(G,G1,G2) :: \int_x G*phi_i^*phi_j ds -\int_x G1*phi_i*d/dx(G2*phi_j) dx
% FuncLap(G,G1,G2) :: \int_x G*phi_i^'*phi_j ds -\int_x G1*phi_i'*d/dx(G2*phi_j) dx
% term1 and term2 denote electric field component 'EComp' 
% eE(xi * df/dp + 1-xi^2 / p * df/dxi)
%   term1 = [formT1,formT2]
%   formT1 = FuncMass :: \int_xi xi*phi_i*phi_j dxi
%   formT2 = FuncGrad :: \int_p phi_i'*phi_j dp
%   term2 = [formT1,formT2]
%   formT1 = FuncMass :: \int_p 1/p*phi_i*phi_j dp
%   formT2 = FuncGrad :: \int_xi (1-xi^2)*phi_i'*phi_j dxi
%
% term3 and term4 denote radiation component 'RadComp'
% -1/p^2 d/dp(stuff*f) + d/dxi(stuff*f)
% 
% term5,6, and 7 denote collision operator component 'CollComp'
% C{F}
%
% Note: All terms are constructed by weak forms
%--------------------------------------------------
Dim = 2;
B = 10000;
tauR=re_tauR(B);
C1 = 1.6e-19 .* E(x,t,params);
C2 = -1/tauR;
ne = 1e19;
Te = 10;
Zeff = 1.4;
phys.c = 2.99792458*10.^8; %speed of light(m/s)
phys.me = 9.10938356*10.^(-31); %mass of electron (kg)
phys.me_keV = 511.000;
phys.e = 1.60217662*10.^(-19); %charge of electron (C)
phys.eo = 8.854187817*10.^(-12); %epsilon0 (C/Vm) or (F/m)
phys.kb = 1.38064852e-23;

%Ecomp
formT1.dim = 2;
formT1.type = 'FuncMass';
formT1.G = @(x)x; %xi
formT2.dim = 1;
formT2.type = 'FuncGrad';
formT2.G = @(x)1; %p
term1 = {formT1,formT2};

formT1.dim = 1;
formT1.type = 'FuncMass';
formT1.G = @(x)1/x; %1/p
formT2.dim = 2;
formT2.type = 'FuncGrad';
formT2.G = @(x)(1-x.^2); %xi
term2 = {formT1,formT2};

%RadComp
formT1.dim = 2;
formT1.type = 'FuncMass';
formT1.G = @(x)(1-x.^2);%xi
formT2.dim = 1;
formT2.type = 'FuncGrad2';
formT2.G = @(x)re_lowerGamma(re_v(x)).*x;%lowerGamma*p
formT2.G1 = @(x)re_lowerGamma(re_v(x)).*x.^3;%lowerGamma*p^3
formT2.G2 = @(x)1/x.^2;%1/p^2
term3 = {formT1,formT2};

formT1.dim = 1;
formT1.type = 'FuncMass';
formT1.G = @(x)1/re_lowerGamma(re_v(x));%1/lowerGamma
formT2.dim = 2;
formT2.type = 'FuncGrad';
formT2.G = @(x)x*(1-x.^2);%xi
term4 = {formT1,formT2};

%CollComp
formT1.dim = 2;
formT1.type = 'FuncMass';
formT1.G = @(x)1;%xi
formT2.dim = 1;
formT2.type = 'FuncGrad2';
formT2.G = @(x)re_C_F(Te,   re_Gamma(ne,re_coulombLog(ne,Te)),   re_Psi(re_x(re_v(x),re_vThermal(Te)))); %C_F
formT2.G1 = @(x)re_C_F(Te,   re_Gamma(ne,re_coulombLog(ne,Te)),   re_Psi(re_x(re_v(x),re_vThermal(Te)))).*x.^2;%p^2*C_F
formT2.G2 = @(x)1/x.^2;%1/p^2
term5 = {formT1,formT2};

formT1.dim = 2;
formT1.type = 'FuncMass';
formT1.G = @(x)1;%xi
formT2.dim = 1;
formT2.type = 'FuncGrad2';
formT2.G = @(x)re_C_A(re_Gamma(ne,re_coulombLog(ne,Te)),   re_v(x),   re_Psi(re_x(re_v(x),re_vThermal(Te))));%C_A
formT2.G1 = @(x)re_C_A(re_Gamma(ne,re_coulombLog(ne,Te)),   re_v(x),   re_Psi(re_x(re_v(x),re_vThermal(Te)))).*x.^2;%C_A*p^2
formT2.G2 = @(x)1/x.^2;%1/p^2
term6 = {formT1,formT2};

formT1.dim = 1;
formT1.type = 'FuncMass';
formT1.G = @(x)re_C_B(Zeff,  re_Gamma(ne,re_coulombLog(ne,Te)),  re_delta(re_vThermal(Te)),  re_v(x),  re_x(re_v(x),re_vThermal(Te)),  re_Psi(re_x(re_v(x),re_vThermal(Te))))/x.^2;%C_B/p^2
formT2.dim = 2;
formT2.type = 'FuncLap';
formT2.G = @(x)(1-x.^2);%xi
formT2.G1 = @(x)(1-x.^2);%xi
formT2.G2 = @(x)1;
term7 = {formT1,formT2};

terms = {term1,term2,term3,term4,term5,term6,term7};

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
function[C_A] = re_C_A(Gamma,v,Psi)
  C_A = Gamma./v.*Psi;
endfunction
function[C_B] = re_C_B(Zeff,Gamma,delta,v,x,Psi)
  C_B = Gamma./(2 .*v).*(Zeff .+erf(x).-Psi.+delta.^4.*x.^2 ./2);
endfunction
function[C_F] = re_C_F(Te,Gamma,Psi)
  C_F = Gamma./Te.*Psi;
endfunction
function [Gamma] = re_Gamma(ne,coulombLog)
  %used to calculate uppercase Gamma from the electron density (ne) and the background electron Temperature (Te)
  Gamma = ne*phys.e.^4*coulombLog/(4*pi*phys.eo.^2);
endfunction
function [coulombLog] = re_coulombLog(ne,Te)
  %used to calculate the coulombLog
  coulombLog = 24 .-log(sqrt(ne/1e6)./(Te*1e3));
endfunction
function [lowerGamma] = re_lowerGamma(v)
  %used to calculate the relativistic gamma if given velocity
  lowerGamma = 1./sqrt(1-rdivide(v.^2, phys.c.^2)); 
endfunction
function [delta] = re_delta(ve)
  %used to calculate delta from the background electron Temperature (Te)
  delta = ve/phys.c;
endfunction
function[ve]=re_vThermal(Te)
  ve = sqrt(2*Te*phys.e/phys.me);
endfunction
function [v] = re_v(rho)
  %used to calculate velocity if given momentum
  v = rho./sqrt((rho./phys.c).^2+phys.me.^2);
endfunction
function [x] = re_x(v,ve)
  x=v/ve;
endfunction
function [Psi] = re_Psi(x)
  %used to calculate Psi
  Psi = 1 ./(2 .*x.^2).*(erf(x) .-x .*(2 ./sqrt(pi) .*exp(-x.^2))); 
endfunction
function[tauR]=re_tauR(B)
  tauR = 6*pi*phys.eo*(phys.me*phys.c).^3/phys.e.^4/B.^2;
endfunction