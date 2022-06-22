function [A] = kronrealspace2DtoDG(pde,opts,deg_out)
%Determines sparse matrix transformation from kronecker 2D realspace
% [phi_1(x),phi_2(x),...,phi_M(x)] \otimes [psi_1(v),psi_2(v),...,psi_N(v)]
% to cell local DG basis where in a cell order is 
% [phi_{l,1}(x),phi_{l,2}(x),...,phi_{l,k+1}(x)] 
%    \otimes
% [psi_{l,1}(v),psi_{l,2}(v),...,psi_{l,k+1}(v)] 
%    for l is the cell index bases on cell indices of (x,v);

% deg_out specifies what order of polynomials you want back.  For example,
% if you only want cell averages, then pass deg_out = 1;

if nargin < 3
    deg_out = opts.deg;
end


%%ONLY extracts constants for now...will add other polynomials later

num_x = 2^pde.dimensions{1}.lev;
num_v = 2^pde.dimensions{2}.lev;
deg = opts.deg;

FG_dof = num_x*num_v*deg^2;
%EYE = speye(FG_dof);

loc_idx = zeros(deg_out^2,1);
for k=0:deg_out-1
    loc_idx(k*deg_out+1:(k+1)*deg_out) = (0:deg_out-1)+deg*num_v*k;
end
    

%AA = sparse(num_x*num_v*deg_out^2,FG_dof);
I = zeros(num_x*num_v*deg_out^2,1);
J = zeros(num_x*num_v*deg_out^2,1);
S = zeros(num_x*num_v*deg_out^2,1);
count = 0;
for i = 1:num_x
    for j=1:num_v
        start = (i-1)*deg^2*num_v+(j-1)*deg;
        %A((i-1)*num_v+1:i*num_v,:) = I(start+1+(0:num_v-1)*deg,:);
        %AA(count+1:count+deg_out^2,:) = EYE(start+1+loc_idx,:);
        %[II,JJ,SS] = find(EYE(start+1+loc_idx,:));
        I(count+1:count+deg_out^2) = count+1:count+deg_out^2;
        J(count+1:count+deg_out^2) = start+1+loc_idx;
        S(count+1:count+deg_out^2) = ones(deg_out^2,1);
        count = count + deg_out^2;
    end
end

A = sparse(I,J,S,num_x*num_v*deg_out^2,FG_dof);


end

