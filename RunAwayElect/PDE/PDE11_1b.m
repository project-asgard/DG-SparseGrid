nuEE = 1;
vT = 1;
delta = 0.042;
Z = 1;
E = 0.0025;
tau = 10^5;
gamma = @(p)sqrt(1+(delta*p).^2);
vx = @(p)1/vT*(p./gamma(p));

Ca = @(p)nuEE*vT^2*(FunPsi(vx(p))./vx(p));

Cb = @(p)1/2*nuEE*vT^2*1./vx(p).*(Z+FunPhi(vx(p))-FunPsi(vx(p))+delta^4*vx(p).^2/2);

Cf = @(p)2*nuEE*vT*FunPsi(vx(p));

Q = 1;

% Exa0 = @(xi,p,t) exp(-2/delta^2*sqrt(1+delta^2*p.^2))/8.5318e-02;
N = 1000;
h = 20/N;
Q = 0;
Fun = @(p)exp(-2/delta^2*sqrt(1+delta^2*p.^2));
for i = 1:N
    x0 = (i-1)*h;
    x1 = i*h;
    [xi,w] = lgwt(20,x0,x1);
    Q = Q+sum(w.*Fun(xi).*xi.^2);
end


% Exa0 = @(xi,p,t) exp(-2/delta^2*sqrt(1+delta^2*p.^2))/(2*Q);
% Exa1 = @(xi,p,t)2/sqrt(pi)*exp(-p.^2);

a = 2;
Exa0 = @(xi,p,t)exp(-(p).^2/a^2)*2/(sqrt(pi)*a^3);

%  Exa0 = @(xi,p,t)(ExactF2(p,t));

fp_bcL = 1; fp_bcR = 0;
fp_bcR = 1;% bR for p
qp_bcL = 0; qp_bcR = 1;
qp_bcR = 0;



fx_bcL = 1; fx_bcR = 1;
qx_bcL = 0; qx_bcR = 0;

function y = FunPsi(x)
psi = @(x)(1./(2*x.^2).*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));
y = psi(x);
ix = find(abs(x)<1e-5);
y(ix) = 0;

end

function y = FunPhi(x)

phi = @(x)erf(x);
y = phi(x);


end

function y = Exa(x)
% y = 3*(1-0.5*sin(pi*x));%exp(-x.^2);
% h0 = 3; h1 = 0.5; h2 = 1; h3 = 0.7; h4 = 3; h6 = 3;
% % ExactF = @(x,t)(h0+h1*x+h2*(3*x.^2-1)/2+h3*(5*x.^3-3*x)/2+...
% %     h4*(35*x.^4-30*x.^2+3)/8+h6*(231*x.^6-315*x.^4+105*x.^2-5)/16);
% t = 0;
y = ...
(h0+h1*x*exp(-1*(1+1)*t)+h2*(3*x.^2-1)/2*exp(-2*(2+1)*t)+h3*(5*x.^3-3*x)/2*exp(-3*(3+1)*t)+...
    h4*(35*x.^4-30*x.^2+3)/8*exp(-4*(4+1)*t)+h6*(231*x.^6-315*x.^4+105*x.^2-5)/16*exp(-6*(6+1)*t));

end

function y = ExactF2(x,t)
y = x-x + (3/5^3)/2;
ix = find(x>5);
y(ix) = 0;
end