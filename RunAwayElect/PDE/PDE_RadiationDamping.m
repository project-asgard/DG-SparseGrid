%This code presents the analytical solution of the Radiation Damping
sigma = 0.1;
f0 = @(x)( exp(-(x-0.36).^2/sigma^2) );
% f0 = @(x)( exp(-(x-0.36).^2/sigma^2)+exp(-(x+0.36).^2/sigma^2) );
phi = @(x,t)( x.*exp(-t)./sqrt(1+(exp(-2*t)-1).*x.^2) );
ExactF = @(x,t)(...
    (phi(x,t).*(1-phi(x,t).^2))./(x.*(1-x.^2)).*f0(phi(x,t)) ...
    );
source = @(x,t)(x-x);


E = 0;
C = 0;
R = 1;