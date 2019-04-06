function flops = kron_mult_cost4( A1,A2,A3,A4)
% flops = kron_mult_cost4( A1,A2,A3,A4)
%
nkron = 4;
rc = zeros(2,nkron);
rc(1,1) = size(A1,1); rc(2,1) = size(A1,2);
rc(1,2) = size(A2,1); rc(2,2) = size(A2,2);
rc(1,3) = size(A3,1); rc(2,3) = size(A3,2);
rc(1,4) = size(A4,1); rc(2,4) = size(A4,2);

nz = zeros(nkron,1);
nz(1) = nnz(A1);
nz(2) = nnz(A2);
nz(3) = nnz(A3);
nz(4) = nnz(A4);

xsizes = prod( rc(2,:) );

flops = kron_mult_cost( rc,nz,xsizes);

