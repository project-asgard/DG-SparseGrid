function [ A ] = kron_realspace_to_DG( lev_vec, deg )
%Determines sparse matrix transformation from kronecker 2D realspace
% [phi_1(x),phi_2(x),...,phi_M(x)] \otimes [psi_1(v),psi_2(v),...,psi_N(v)]
% to cell local DG basis where in a cell order is 
% [phi_{l,1}(x),phi_{l,2}(x),...,phi_{l,k+1}(x)] 
%    \otimes
% [psi_{l,1}(v),psi_{l,2}(v),...,psi_{l,k+1}(v)] 
%    for l is the cell index bases on cell indices of (x,v);

num_dims = numel(lev_vec);

assert( num_dims <= 2, 'kron_realspace_to_DG only implemented for num_dims<=2' )

if( num_dims == 1 )

    num_x = 2^lev_vec(1);

    A = sparse( eye(num_x*deg) );

else

    num_x = 2^lev_vec(1);
    num_v = 2^lev_vec(2);

    FG_dof = num_x*num_v*deg^2;

    loc_idx = zeros(deg^2,1);
    for k=0:deg-1
        loc_idx(k*deg+1:(k+1)*deg) = (0:deg-1)+deg*num_v*k;
    end
    
    I = zeros(num_x*num_v*deg^2,1);
    J = zeros(num_x*num_v*deg^2,1);
    S = zeros(num_x*num_v*deg^2,1);
    count = 0;
    for i = 1:num_x
        for j=1:num_v
            start = (i-1)*deg^2*num_v+(j-1)*deg;
            I(count+1:count+deg^2) = count+1:count+deg^2;
            J(count+1:count+deg^2) = start+1+loc_idx;
            S(count+1:count+deg^2) = ones(deg^2,1);
            count = count + deg^2;
        end
    end

    A = sparse(I,J,S,num_x*num_v*deg^2,FG_dof);

end

end
