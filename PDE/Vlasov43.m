function pde=Vlasov43
% Example 4.3for Vlasov Equation
% Two-stream instability

% parameters

p.k_0 = 2.0/13.0;
p.A = 0.05;
p.v_t = 0.3;
p.u = 0.99;

p.Lmin = 0;
p.Lmax = 13*pi;
p.Vmin = -5;
p.Vmax = +5;

p.TEND = 1;

pde.Fx_0 = @Fx_0;
pde.Fv_0 = @Fv_0;
pde.Fxv_0 = @Fxv_0;
pde.Ex = @Ex;
pde.Et = @Et;
pde.E = @E;
pde.p = p;

pde.params = p;

pde.solvePoisson = 1;
pde.applySpecifiedE = 0;
pde.implicit = 0;
pde.checkAnalytic = 0;

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
A = p.A;
k_0 = p.k_0;
f=(1+A*cos(k_0*x));
end
function f=Fv_0(v,p)
% Initial condition for v variable
u = p.u;
v_t = p.v_t;
f=1/(2*v_t*sqrt(2*pi)) * ( exp(-abs(u+v).^2/(2*v_t^2)) + exp(-abs(u-v).^2/(2*v_t^2)) );
end
function f=Fxv_0(x,v,p)
f=Fv_0(v).*Fx_0(x);
end

function f=Ex(x,p)
f=x.*0;
end
function f=Et(x,p)
f=x.*0;
end
function f=E(x,t,p)
f=x.*0;
end

% source term--fully seperable functions
% source = source1+source2+source3

% source term 1
function f = source1t(t)
f = t.*0;
end
function f = source1x(x,p)
f = x.*0;
end
function f = source1v(v,p)
f = v.*0;
end
function f = source1(x,v,t)
f = source1x(x).*source1v(v).*source1t(t);
end

% source term 2
function f = source2t(t)
f = t.*0;
end
function f = source2x(x,p)
f = x.*0;
end
function f = source2v(v,p)
f = v.*0;
end
function f = source2(x,v,t)
f = source2x(x).*source2v(v).*source2t(t);
end

% source term 3
function f = source3t(t)
f = t.*0;
end
function f = source3x(x,p)
f = x.*0;
end
function f = source3v(v,p)
f = v.*0;
end
function f = source3(x,v,t)
f = source3x(x).*source3v(v).*source3t(t);
end


% Exact F
function f=ExactFt(t)
f=t.*0;
end
function f=ExactFx(x,p)
f = x.*0;
end
function f=ExactFv(v,p)
f = v.*0;
end
function f=ExactF(x,v,t)
f = ExactFx(x).*ExactFv(v).*ExactFt(t);
end

