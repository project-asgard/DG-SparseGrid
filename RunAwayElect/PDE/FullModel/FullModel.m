nuEE = 1;
vT = 1;
delta = 0.042;
Z = 1;
E = 0.0025;
tau = 10^5;
gamma = @(p)sqrt(1+(delta*p).^2);
vx = @(p)1/vT*(p./gamma(p));

% phi = @(x)erf(x);
% psi = @(x)1./(2*x.^2).*(phi(x)-2*x/sqrt(pi).*exp(-x.^2));

Ca = @(p)nuEE*vT^2*(FunPhi(vx(p))./vx(p));%@(p)nuEE*vT^2*(phi(vx(p))./vx(p));

Cb = @(p)1/2*nuEE*vT^2*1./vx(p).*(Z+FunPhi(vx(p))-FunPsi(vx(p))+delta^4*vx(p).^2/2);

Cf = @(p)2*nuEE*vT*FunPsi(vx(p));

Exa0 = @(p,xi,t)2/sqrt(pi)*exp(-p.^2);

function y = FunPsi(x)
psi = @(x)(1./(2*x.^2).*(erf(x)-2*x/sqrt(pi).*exp(-x.^2)));
y = psi(x);
ix = find(abs(x)<1e-5);
y(ix) = 0;

end

function y = FunPhi(x)

phi = @(x)erf(x);
y = phi(x);
ix = find(abs(x)<1e-5);
y(ix) = 0;
% if abs(x)<1e-5
%     y = 0;
% else
%     y = x.*psi(x);
% end
end