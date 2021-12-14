function [pterm] = pterm_coeff_matrix(deg,t,dim,pterm,params, ...
    FMWT_blocks, coeff_level)

% Get the RHS matrix (L) for this pterm

[L,L_not_rotated] = coeff_matrix(deg,t,dim,pterm,params,FMWT_blocks,coeff_level); % note that dV is applied in here

% Get the LHS mass matrix (M) for this pterm
% this mass matrix needs to have the volume jacobian incorporated instead
% of the possibly different surface jacobian.
lhs_mass_pterm = MASS(pterm.LHS_mass_g,'','',dim.moment_dV);
[M,M_not_rotated] = coeff_matrix(deg,t,dim,lhs_mass_pterm,params,FMWT_blocks,coeff_level); 

% Move M to the RHS

pterm.mat = inv(M) * L;
pterm.mat_unrotated = inv(M_not_rotated) * L_not_rotated;


% Store M for use in application to the BC

pterm.LHS_mass_mat = M;
pterm.LHS_mass_mat_unrotated = M_not_rotated;

end