function [term_1D] = sd_term_coeff_matrix(deg,t,dim,term_1D,params,FMWT_blocks,coeff_level)

if ~exist('lev','var') || isempty(coeff_level)
    coeff_level = dim.lev;
end

% Get the matrix for each pterm

for i=1:numel(term_1D.pterms)
    term_1D.pterms{i} = pterm_coeff_matrix(deg,t,dim,term_1D.pterms{i},params,FMWT_blocks,coeff_level);
end

% Chain together the partial terms

mat = eye(size(term_1D.pterms{1}.mat));
mat_unrotated = mat;

for i=1:numel(term_1D.pterms)
    mat = mat * term_1D.pterms{i}.mat;
    mat_unrotated = mat_unrotated * term_1D.pterms{i}.mat_unrotated;
end

term_1D.mat = mat;
term_1D.mat_unrotated = mat_unrotated;

end


