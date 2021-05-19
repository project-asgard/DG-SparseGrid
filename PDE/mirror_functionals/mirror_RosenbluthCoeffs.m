function vec = mirror_RosenbluthCoeffs(f_prev, lIndex, uval)

%evaluation of A_l (u) and B_l(u), the velocity-space coefficients for the
%Rosenbluth potentials

% A_l (u) = 4pi/(2l + 1) [ 1/(2l+3) (E_l (U_l/gamma) * (u/gamma)^(-l + 1) + M_l
% (U_l/gamma) * (u/gamma)^(l+2) ) - 1/(2l-1) (N_l (U_l/gamma) *
% (u/gamma)^(1-l) + R_l (U_l/gamma) * (u/gamma)^l ) ]

% B_l (u) = 4pi/(2l + 1) [ 1/(2l+3) (E_l (U_l) * (u/gamma)^(-l + 1) + M_l
% (U_l) * (u/gamma)^(l+2) ) - 1/(2l-1) (N_l (U_l) *
% (u/gamma)^(1-l) + R_l (U_l) * (u/gamma)^l ) ]


%speed of light

c = 3*10^8; %m/s

%relativistic correction

gamma = @(u) sqrt(1 + u.^2./c^2);

%argument for functionals in A is U/gamma

U_coeff_A = @(u) mirror_Ucoeff(f_prev, lIndex, u)./gamma(u);

%argument for functions in B is U

U_coeff_B = @(u) mirror_Ucoeff(f_prev, lIndex, u);

%evaluating functionals

funcs_A = mirror_functional(U_coeff_A, uval, lIndex);
funcs_B = mirror_functional(U_coeff_B, uval, lIndex);

%writing out formulae for A and B

A = (4*pi/(2*lIndex + 1))*( (1/(2*lIndex+3))*(funcs_A(4) * (uval/gamma(uval))^(-lIndex + 1) + funcs_A(1) * (uval/gamma(uval))^(lIndex+2) ) - (1/(2*lIndex-1))*(funcs_A(2)*(uval/gamma(uval))^(1-lIndex) + funcs_A(3)*(uval/gamma(uval))^lIndex ) );
B = (4*pi/(2*lIndex + 1))*( (1/(2*lIndex+3))*(funcs_B(4) * (uval/gamma(uval))^(-lIndex + 1) + funcs_B(1) * (uval/gamma(uval))^(lIndex+2) ) - (1/(2*lIndex-1))*(funcs_B(2)*(uval/gamma(uval))^(1-lIndex) + funcs_B(3)*(uval/gamma(uval))^lIndex ) );

vec = [A,B];

end