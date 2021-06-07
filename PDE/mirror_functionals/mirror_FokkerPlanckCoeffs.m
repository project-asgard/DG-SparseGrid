function vec = mirror_FokkerPlanckCoeffs(f_b,uVal,lIndex,z)

%evaluation of local Fokker-Planck coefficients A through F corresponding to
%equations 31-36 in Franz's dissertation (Franz, 1993) for a species a
%colliding with background species b

params = mirror_parameters();

val_z = cos(z);
%getting Rosenbluth coefficients for background species at lIndex

[A_b_l, B_b_l] = mirror_RosenbluhCoeffs(f_b, uVal, lIndex);

%argument for functionals in A is U/gamma

U_coeff_A = @(u) mirror_Ucoeff(f_prev, u, lIndex)./gamma(u);

%argument for functions in B is U

U_coeff_B = @(u) mirror_Ucoeff(f_prev, u, lIndex);

%getting M,N,R,E functionals for A and B coefficients

funcs_A = mirror_functional(U_coeff_A, uVal, lIndex);
funcs_B = mirror_functional(U_coeff_B, uVal, lIndex);

%getting delta, delta^2, and delta^3 of A and B Rosenbluth Coefficients

delta_A = delta_func(funcs_A,uVal,lIndex,params);
delta2_A = delta2_func(funcs_A,uVal,lIndex,params);
delta3_A = delta3_func(funcs_A,uVal,lIndex,params);

delta_B = delta_func(funcs_B,uVal,lIndex,params);
delta2_B = delta2_func(funcs_B,uVal,lIndex,params);
delta3_B = delta3_func(funcs_B,uVal,lIndex,params);


% A_c = -0.5*gamma^2*sum_{b,l}[(Z_b/Z_a)^2*(m_a/m_b)*ln_delt(2*l*(l+1)*(gamma/u)*A^b_l
% - (2 + l*(l+1))*(delta_A^b_l) + (u/gamma^3)*2*(delta2_A^b_l +
% (u/gamma)*delta3_A^b_l)]*P_l (cos (theta))
 
%individual element for background species with distribution f_b at order lIndex
A_c = @(u) -0.5*params.gamma(u).^2*[(params.b.Z/params.a.Z)^2*(params.a.m/params.b.m)*paramas.ln_delt.*(2*lIndex*(lIndex+1)*(params.gamma(u)/u).*A_b_l- (2 + lIndex*(lIndex+1))*delta_A + (uVal/params.gamma(uVal)^3)*(2*delta2_A +(uVal/params.gamma(uVal))*delta3_A) ]*(legendre(lIndex,val_z));

% B_c = 0.5*u^2*sum_{b,l}[(Z_b/Z_a)^2*ln_delt*(delta2_B^b_l]*P_l (cos (theta))

B_c = 0.5*uVal^2*[(params.b.Z/params.a.Z)^2*params.ln_delt*delta2_B]*legendre(lIndex,val_z);

% C_c = 0.5*gamma*sum_{b,l}[(Z_b/Z_a)^2*ln_delt*(delta_B^b_l) -
% (gamma/u)*B^b_l]*(dP_l (cos (theta))/dtheta)

C_c = 0.5*params.gamma(uVal)^2*[(params.b.Z/params.a.Z)^2*params.ln_delt*delta_B*(params.gamma(uVal)/uVal)*B_b_l]*(-lIndex-1)*(val_z*legendre(lIndex,val_z) - legendre(lIndex+1,val_z))/(val_z^2 -1);

% D_c = -0.5*gamma^2*sum_{b,l}[(Z_b/Z_a)^2*(m_a/m_b)*ln_delt*(delta2_A^b_l/gamma^3) - 
% (2*delta_A^b_l/u) - l*(l+1)*(gamma/u^2)*A^b_l]*sin (theta) dP_l (cos (theta))/dtheta

D_c = -0.5*params.gamma(uVal)^2*[(params.b.Z/params.a.Z)^2*(params.a.m/params.b.m)*params.ln_delt*(delta2_A/params.gamma(uVal)^3) - (2*delta_A - lIndex*(lIndex+1)*(params.gamma(uVal)/uVal^2)*A_b_l]*sin(z)*(-lIndex-1)*(val_z*legendre(lIndex,val_z) - legendre(lIndex+1,val_z))/(val_z^2 -1);


% E_c = 0.5*gamma*sum_{b,l}[(Z_b/Z_a)^2*ln_delt*(delta_B^b_l) -
% (gamma/u)*B^b_l]*sin(theta)*(dP_l (cos (theta))/dtheta)

E_c = 0.5*params.gamma(uVal)*[(params.b.Z/params.b.Z)^2*params.ln_delt*(delta_B) -(params.gamma(uVal)/uVal)*B_b_l]*sin(z)*(-lIndex-1)*(val_z*legendre(lIndex,val_z) - legendre(lIndex+1,val_z))/(val_z^2 -1);

% F_c = 0.5*sum_{b,l}[(Z_b/Z_a)^2*ln_delt*(delta_B^b_l)*(P_l (cos (theta))
% + (gamma^2/u^2)*B^b_l*(d^2 P_l (cos (theta))/dtheta^2)]*sin(theta)

F_c = 0.5*[(params.b.Z/params.a.Z)^2*params.ln_delt*(delta_B)*(legendre(lIndex,val_z)) + (params.gamma(uVal)^2/uVal^2)*B_b_l*(lIndex+1)*((lIndex*(val_z^2 -1) + 2*val_z^2)*legendre(lIndex,val_z) - 2*val_z*legendre(lIndex+1,val_z))/(val_z^2 -1)^2]*sin(z);

%delta functions from Franz thesis

    function val = delta_func(funcs,uVal,lIndex, params)
        
        val = (4*pi/(2*lIndex + 1))*( (1/(2*lIndex+3))*(funcs(1) *(lIndex+2)*(uVal/params.gamma(uVal))^(lIndex + 1) - funcs(4) * (lIndex+1)*(uVal/params.gamma(uVal))^(-lIndex-2) ) - (1/(2*lIndex-1))*(funcs(3)*lIndex*(uVal/params.gamma(uVal))^(lIndex -1) - funcs(2)*(lIndex-1)*(uVal/params.gamma(uVal))^(-lIndex) ) );
        
    end

    function val1 = delta2_func(funcs,uVal,lIndex, params)
        
        val1 = (4*pi/(2*lIndex + 1))*( ((lIndex+1)*(lIndex+2)/(2*lIndex+3))*(funcs(1)*(uVal/params.gamma(uVal))^lIndex + funcs(4)*(uVal/params.gamma(uVal))^(-lIndex-3) ) - (lIndex*(lIndex-1)/(2*lIndex-1))*(funcs(3)*(uVal/params.gamma(uVal))^(lIndex -2) - funcs(2)*(uVal/params.gamma(uVal))^(-lIndex-1) ) );
        
    end

    function val2 = delta3_func(funcs,uVal,lIndex, params)
        
        val2 = (4*pi/(2*lIndex + 1))*( ((lIndex+1)*(lIndex+2)/(2*lIndex+3))*(funcs(1)*lIndex*(uVal/params.gamma(uVal))^(lIndex-1) - funcs(4) * (lIndex+3)*(uVal/params.gamma(uVal))^(-lIndex-4) ) - (lIndex*(lIndex-1)/(2*lIndex-1))*(funcs(3)*(lIndex -2)*(uVal/params.gamma(uVal))^(lIndex -3) - funcs(2)*(lIndex+1)*(uVal/params.gamma(uVal))^(-lIndex-2) ) );
        
    end

vec = [A_c,B_c,C_c,D_c,E_c,F_c];


end