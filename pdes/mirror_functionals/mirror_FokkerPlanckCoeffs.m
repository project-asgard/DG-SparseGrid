function vec = mirror_FokkerPlanckCoeffs(f_b,uVal,lIndex,z,params,b_species)

%evaluation of local Fokker-Planck coefficients A through F corresponding to
%equations 31-36 in Franz's dissertation (Franz, 1993) for a species a
%colliding with background species b
global abs_tol; 
global rel_tol;
abs_tol = 1e-8;
rel_tol = 1e-8;

params.b = b_species;

val_z = cos(z);
legendre_val = legendre(lIndex,val_z);
legendre_val = legendre_val(lIndex+1);
legendre_val1 = legendre(lIndex+1,val_z);
legendre_val1 = legendre_val1(lIndex+2);
%getting Rosenbluth coefficients for background species at lIndex

coeffs_R = mirror_RosenbluthCoeffs(f_b, uVal, lIndex);

%argument for functionals in A is U/gamma

U_coeff_A = @(u) mirror_Ucoeff(f_b, u, lIndex)./params.gamma(u);

%argument for functions in B is U

U_coeff_B = @(u) mirror_Ucoeff(f_b, u, lIndex);

%getting M,N,R,E functionals for A and B coefficients

funcs_A = mirror_functional(U_coeff_A, uVal, lIndex);
funcs_B = mirror_functional(U_coeff_B, uVal, lIndex);

%getting delta, delta^2, and delta^3 of A and B Rosenbluth Coefficients

U_coeff_A(uVal)
U_coeff_B(uVal)

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
A_c = -0.5*params.gamma(uVal).^2*((params.b.Z/params.a.Z)^2*(params.a.m/params.b.m)*params.ln_delt* ... 
    (2*lIndex*(lIndex+1)*(params.gamma(uVal)/uVal).*coeffs_R(1) - (2 + lIndex*(lIndex+1))*delta_A ...
    + (uVal/params.gamma(uVal)^3)*(2*delta2_A +(uVal/params.gamma(uVal))*delta3_A) ))*legendre_val;

% B_c = 0.5*u^2*sum_{b,l}[(Z_b/Z_a)^2*ln_delt*(delta2_B^b_l]*P_l (cos (theta))

B_c = 0.5*uVal^2*((params.b.Z/params.a.Z)^2*params.ln_delt*delta2_B)*legendre_val;

% C_c = 0.5*gamma*sum_{b,l}[(Z_b/Z_a)^2*ln_delt*(delta_B^b_l) -
% (gamma/u)*B^b_l]*(dP_l (cos (theta))/dtheta)

C_c = 0.5*params.gamma(uVal)*((params.b.Z/params.a.Z)^2*params.ln_delt*delta_B - (params.gamma(uVal)/uVal)*coeffs_R(2))...
    *legendre_derivative(lIndex,val_z);

% D_c = -0.5*gamma^2*sum_{b,l}[(Z_b/Z_a)^2*(m_a/m_b)*ln_delt*(delta2_A^b_l/gamma^3) - 
% (2*delta_A^b_l/u) - l*(l+1)*(gamma/u^2)*A^b_l]*sin (theta) dP_l (cos
% (theta))/dtheta (-sin(theta))

D_c = -0.5*params.gamma(uVal)^2*((params.b.Z/params.a.Z)^2*(params.a.m/params.b.m)*params.ln_delt*...
    (delta2_A/params.gamma(uVal)^3) - (2*delta_A - lIndex*(lIndex+1)*(params.gamma(uVal)/uVal^2)*coeffs_R(1)))*... 
    legendre_derivative(lIndex,val_z);


% E_c = 0.5*gamma*sum_{b,l}[(Z_b/Z_a)^2*ln_delt*(delta_B^b_l) -
% (gamma/u)*B^b_l]*sin(theta)*(dP_l (cos (theta))/dtheta)

E_c = 0.5*params.gamma(uVal)*((params.b.Z/params.b.Z)^2*params.ln_delt*(delta_B) -(params.gamma(uVal)/uVal)*coeffs_R(2))...
*(sin(z))*(-lIndex-1)*legendre_derivative(lIndex,val_z);

% F_c = ln_delt*sum_{b,l}*0.5*(Z_b/Z_a)^2[(gamma/u)*(delta_B^b_l)*(P_l (cos (theta))
% + (gamma^2/u^2)*B^b_l*(d^2 P_l (cos (theta))/dtheta^2)]*sin(theta)

F_c = 0.5*(params.b.Z/params.a.Z)^2*params.ln_delt*((params.gamma(uVal)/uVal)*(delta_B)*legendre_val ...
+ (params.gamma(uVal)^2/uVal^2)*coeffs_R(2)*legendre_derivative2(lIndex,val_z))*(sin(z));

%derivatives of legendre polynomials

    function val = legendre_derivative(l,val_z)
        if (l == 0)
            val = 0;
        else
            legendre_l = legendre(l,val_z);
            legendre_l = legendre_l(l+1);
            legendre_l1= legendre(l-1,val_z);
            legendre_l1 = legendre_l1(l);
            val = (val_z*legendre_l - legendre_l1)/(val_z^2 -1);
        end
        
    end

    function val = legendre_derivative2(l,val_z)
        if (l == 0)
            val = 0;
        else
            %getting first derivative values
            legendre_d1_l = legendre_derivative(l,val_z);
            legendre_d1_l1 = legendre_derivative(l-1,val_z);
            legendre_l = legendre(l,val_z);
            legendre_l = legendre_l(l+1);
            val = -2*val_z*legendre_d1_l/(val_z^2 -1) + (l/(val_z^2 -1))*(legendre_l + val_z*legendre_d1_l ... 
                - legendre_d1_l1);
        end
        
    end

%delta functions from Franz thesis

    function val = delta_func(funcs,uVal,lIndex, params)
        
        val = (4*pi/(2*lIndex + 1))*( (1/(2*lIndex+3))*(funcs(1) *(lIndex+2)*(uVal/params.gamma(uVal))^(lIndex + 1) ...
            - funcs(4) * (lIndex+1)*(uVal/params.gamma(uVal))^(-lIndex-2) ) - (1/(2*lIndex-1))...
            *(funcs(3)*lIndex*(uVal/params.gamma(uVal))^(lIndex -1) - funcs(2)*(lIndex-1)...
            *(uVal/params.gamma(uVal))^(-lIndex) ) );
        
    end

    function val1 = delta2_func(funcs,uVal,lIndex, params)
        
        val1 = (4*pi/(2*lIndex + 1))*( ((lIndex+1)*(lIndex+2)/(2*lIndex+3))*(funcs(1)...
            *(uVal/params.gamma(uVal))^lIndex + funcs(4)*(uVal/params.gamma(uVal))^(-lIndex-3) ) ...
            - (lIndex*(lIndex-1)/(2*lIndex-1))*(funcs(3)*(uVal/params.gamma(uVal))^(lIndex -2) - ...
            funcs(2)*(uVal/params.gamma(uVal))^(-lIndex-1) ) );
        
    end

    function val2 = delta3_func(funcs,uVal,lIndex, params)
        
        val2 = (4*pi/(2*lIndex + 1))*( ((lIndex+1)*(lIndex+2)/(2*lIndex+3))*(funcs(1)*lIndex ...
            *(uVal/params.gamma(uVal))^(lIndex-1) - funcs(4) * (lIndex+3)*(uVal/params.gamma(uVal))^(-lIndex-4) ) ... 
            - (lIndex*(lIndex-1)/(2*lIndex-1))*(funcs(3)*(lIndex -2)*(uVal/params.gamma(uVal))^(lIndex -3) ... 
            - funcs(2)*(lIndex+1)*(uVal/params.gamma(uVal))^(-lIndex-2) ) );
        
    end

vec = [A_c,B_c,C_c,D_c,E_c,F_c];


end