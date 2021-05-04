function ans = mirror_Rfunctional(wfunc,umin,lindex)

% evaluation of R_l (Franz, 1987) using input function w(u) for value of
% index value lindex
% output is a 1D function in u-space

% R_l(w) = int^{u_max}_u w(u') (u'/gamma')^{3-l} gamma'^2 du'

%speed of light

c = 3*10^8; %m/s

%setting up integrand matrix

gamma = @(u) sqrt(1 + u.^2./c^2); %relativistic correction

R_int = @(u) gamma(u).^2.*(u./gamma(u)).^(3 - lindex);

%M_int = reshape(M_int,size(wfunc));
%gamma = reshape(gamma, size(wfunc));
integrand = @(u) wfunc(u).*R_int(u);
%integrating
ans = integral(integrand, umin, Inf);

end