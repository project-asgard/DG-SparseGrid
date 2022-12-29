function [pterm] = pterm_coeff_matrix(deg,t,dim,pterm,params, ...
    FMWT_blocks, coeff_level)

% Get the RHS matrix (L) for this pterm

[L,L_not_rotated] = coeff_matrix(deg,t,dim,pterm,params,FMWT_blocks,coeff_level); % note that dV is applied in here

% Get the LHS mass matrix (M) for this pterm
% this mass matrix needs to have the volume jacobian incorporated instead
% of the possibly different surface jacobian.
if isempty(pterm.LHS_mass_mat)
    lhs_mass_pterm = MASS(pterm.LHS_mass_g,'','',dim.moment_dV);
    [M,M_not_rotated] = coeff_matrix(deg,t,dim,lhs_mass_pterm,params,FMWT_blocks,coeff_level); 
    pterm.LHS_mass_mat = full(M);
    pterm.LHS_mass_mat_unrotated = full(M_not_rotated);
    identity_ref = norm(pterm.LHS_mass_mat-eye(size(M)),'fro');
    if identity_ref < 1e-11 %%Checking for identity mass matrix
        %fprintf('Replacing mass mat with identity.  Error = %e\n',identity_ref);
        pterm.LHS_mass_mat = speye(size(M));
        pterm.LHS_mass_mat_unrotated = speye(size(M));
    end
end

% Move M to the RHS

pterm.mat = pterm.LHS_mass_mat\L;
pterm.mat_unrotated = pterm.LHS_mass_mat_unrotated\L_not_rotated;


% Store M for use in application to the BC

end