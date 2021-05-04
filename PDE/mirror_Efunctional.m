function ans = mirror_Efunctional(wfunc,umax,lindex)

% evaluation of E_l (Franz, 1987) using input function w(u) for value of
% index value lindex
% output is a 1D function in u-space

% E_l(w) = int^{u_max}_u w(u') (u'/gamma')^{1-l} gamma'^2 du'

%speed of light

c = 3*10^8; %m/s

%setting up integrand matrix

gamma = @(u) sqrt(1 + u.^2./c^2); %relativistic correction

M_int = @(u) gamma(u).^2.*(u./gamma(u)).^(4 + lindex);

%M_int = reshape(M_int,size(wfunc));
%gamma = reshape(gamma, size(wfunc));
integrand = @(u) wfunc(u).*M_int(u);
%integrating
ans = integral(integrand, 0, umax);

end