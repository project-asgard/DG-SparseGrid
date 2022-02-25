function [x] = multigrid(opts,pde,A_data,k,lam,b,x,imex_flag)
%Multigrid preconditioner for bicgstabl

if nargin < 8
    x = zeros(size(b));
end

x = twogrid(opts,pde,A_data,k,lam,b,x);

end

function [x] = twogrid(opts,pde,A_data,k,lam,b,x,imex_flag)

    if numel(b) == 1
        x = b/fast_2d_matrix_apply(opts,pde,A_data,1,imex_flag);
    else
        
    end

end
