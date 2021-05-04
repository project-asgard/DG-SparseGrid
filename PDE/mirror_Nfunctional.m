function ans = mirror_Nfunctional(wfunc,umax,lindex)

% evaluation of N_l (Franz, 1987) using input function w(u) for value of
% index value lindex
% output is a 1D function in u-space

% N_l(w) = int^{u_max}_u w(u') (u'/gamma')^{2+l} gamma'^2 du'

%speed of light

c = 3*10^8; %m/s

%setting up integrand matrix

gamma = @(u) sqrt(1 + u.^2./c^2); %relativistic correction

N_int = @(u) gamma(u).^2.*(u./gamma(u)).^(2 + lindex);

%M_int = reshape(M_int,size(wfunc));
%gamma = reshape(gamma, size(wfunc));
integrand = @(u) wfunc(u).*N_int(u);
%integrating
ans = integral(integrand, 0, umax);

end