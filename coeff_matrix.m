%% Construct 1D coefficient matrices
% This routine returns a 2D array representing an operator coefficient
% matrix for a single dimension (1D). Each term in a PDE requires D many coefficient
% matricies. These operators can only use the supported types below.


function [term_1D] = coeff_matrix(deg,t,dim,term_1D,params,FMWT_blocks,lev)

if ~exist('lev','var') || isempty(lev)
    lev = dim.lev;
end


%%
% Get the matrix for each partial terms

for i=1:numel(term_1D.pterms)
    
    pterm = term_1D.pterms{i};
    [mat,mat1,mat2,mat0] = coeff_matrix_mass_or_grad(deg,t,dim,pterm,params,FMWT_blocks,lev);
    
    term_1D.pterms{i}.mat = mat;
    term_1D.pterms{i}.mat_unrotated = mat0;
    
end

%%
% For non-unitary jacobian systems we need to invert the jacobian mass
% matrix for each of the extra (more than 1) pterms

jacobian_pterm = MASS(dim.jacobian);
mass_J = coeff_matrix_mass_or_grad(deg,t,dim,jacobian_pterm,params,FMWT_blocks,lev);

%%
% Chain together the partial terms
dof = deg * 2^dim.lev;
[m,n] = size(mat);
assert(m==n);
assert(m >= dof);
mat = eye(dof);
mat_unrotated = eye(dof);

for i=1:numel(term_1D.pterms)
    if i>1
        %disp('mult by J^-1');
        mat = mat * inv(mass_J);
    end
    %disp('mult by mat')
    mat = mat * term_1D.pterms{i}.mat(1:dof, 1:dof);
    [m,n] = size(mat); % these three lines are just here to test if I can remove all the (1:dof,1:dof) cruft
    assert(m==dof);
    assert(n==dof)
    mat_unrotated = mat_unrotated * term_1D.pterms{i}.mat_unrotated(1:dof, 1:dof);
end

term_1D.mat = mat;
term_1D.mat_unrotated = mat_unrotated;

end


