function vec = mirror_RosenbluthCoeffs(f_prev, uVal, lIndex)

%evaluation of A_l (u) and B_l(u), the velocity-space coefficients for the
%Rosenbluth potentials

% A_l (u) = 4pi/(2l + 1) [ 1/(2l+3) (E_l (U_l/gamma) * (u/gamma)^(-l + 1) + M_l
% (U_l/gamma) * (u/gamma)^(l+2) ) - 1/(2l-1) (N_l (U_l/gamma) *
% (u/gamma)^(1-l) + R_l (U_l/gamma) * (u/gamma)^l ) ]

% B_l (u) = 4pi/(2l + 1) [ 1/(2l+3) (E_l (U_l) * (u/gamma)^(-l + 1) + M_l
% (U_l) * (u/gamma)^(l+2) ) - 1/(2l-1) (N_l (U_l) *
% (u/gamma)^(1-l) + R_l (U_l) * (u/gamma)^l ) ]


%speed of light

c = 3*10^10; %cm/s

%relativistic correction

gamma = @(u) sqrt(1 + u.^2./c^2);

%argument for functionals in A is U/gamma

U_coeff_A = @(u) mirror_Ucoeff(f_prev, u, lIndex)./gamma(u);

%argument for functions in B is U

U_coeff_B = @(u) mirror_Ucoeff(f_prev, u, lIndex);

%evaluating functionals

funcs_A = mirror_functional(U_coeff_A, uVal, lIndex);
funcs_B = mirror_functional(U_coeff_B, uVal, lIndex);

%writing out formulae for A and B

A = (4*pi/(2*lIndex + 1))*( (1/(2*lIndex+3))*(funcs_A(4) * (uVal/gamma(uVal))^(-lIndex - 1) + funcs_A(1) * (uVal/gamma(uVal))^(lIndex+2) ) - (1/(2*lIndex-1))*(funcs_A(2)*(uVal/gamma(uVal))^(1-lIndex) + funcs_A(3)*(uVal/gamma(uVal))^lIndex ) );
B = (4*pi/(2*lIndex + 1))*( (1/(2*lIndex+3))*(funcs_B(4) * (uVal/gamma(uVal))^(-lIndex - 1) + funcs_B(1) * (uVal/gamma(uVal))^(lIndex+2) ) - (1/(2*lIndex-1))*(funcs_B(2)*(uVal/gamma(uVal))^(1-lIndex) + funcs_B(3)*(uVal/gamma(uVal))^lIndex ) );

vec = [A,B];

end