function vec = mirror_functional(wfunc,uLimit,lIndex)

% evaluation of M_l (Franz, 1987) using input function w(u) for value of
% index value lindex
% output is a 1D function in u-space

% M_l(w) = int^{u_max}_u w(u') (u'/gamma')^{1-l} gamma'^2 du'

%speed of light

c = 3*10^8; %m/s

%setting up integrand matrix

gamma = @(u) sqrt(1 + u.^2./c^2); %relativistic correction

M_int = @(u) gamma(u).^2.*(u./gamma(u)).^(1 - lIndex);
N_int = @(u) gamma(u).^2.*(u./gamma(u)).^(2 + lIndex);
R_int = @(u) gamma(u).^2.*(u./gamma(u)).^(3 - lIndex);
E_int = @(u) gamma(u).^2.*(u./gamma(u)).^(4 + lIndex);


%M_int = reshape(M_int,size(wfunc));
%gamma = reshape(gamma, size(wfunc));
integrand_M = @(u) wfunc(u).*M_int(u);
integrand_N = @(u) wfunc(u).*N_int(u);
integrand_R = @(u) wfunc(u).*R_int(u);
integrand_E = @(u) wfunc(u).*E_int(u);

%integrating
M = integral(integrand_M, uLimit, Inf);
N = integral(integrand_N, 0, uLimit);
R = integral(integrand_R, uLimit, Inf);
E = integral(integrand_E, 0, uLimit);

vec = [M,N,R,E];

end