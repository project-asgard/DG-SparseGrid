function flops = kron_mult_cost2( A1,A2)
% flops = kron_mult_cost2( A1,A2)
%
nkron = 2;
rc = zeros(2,nkron);
rc(1,1) = size(A1,1); rc(2,1) = size(A1,2);
rc(1,2) = size(A2,1); rc(2,2) = size(A2,2);

nz = zeros(nkron,1);
nz(1) = nnz(A1);
nz(2) = nnz(A2);

xsizes = prod( rc(2,:) );

flops = kron_mult_cost( rc,nz,xsizes);

