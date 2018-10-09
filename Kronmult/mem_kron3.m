function [total_mem_use, mem_use] = mem_kron3( nr1, nc1, nr2,nc2, nr3,nc3, nvec )
%
% [total_mem_use, mem_use] = mem_kron3( nr1, nc1, nr2,nc2, nr3,nc3, nvec )
% 
% 
% estimate total amount of temporary space needed
%
idebug = 1;
%  Y =  { kron(A2,A3) * X } * transpose(A1)
%
%  X is (nc3 by nc2 by nc1) * nvec
%  computed as
%  step 1: Ytmp = kron2( A2,A3, X )
%  step 2: Y = Ytmp * transpose(A1)
%
%  note Ytmp is  (nr3 * nr2) by nvec2
%
nvec2 = nc1 * nvec;
mem_use = (nr3 * nr2) * nvec2;

[total_mem_use2, mem_use2] = mem_kron2( nr2,nc2, nr3,nc3, nvec2 );

% ------------------------------------------
% no need for extra temporary array in kron1
% since it is just calling GEMM
% ------------------------------------------
total_mem_use = total_mem_use2 + mem_use;

if (idebug >= 1),
  disp(sprintf('mem_kron3: mem_use=%g, total_mem_use=%g, total_mem_use2=%g, mem_use2=%g', ...
                           mem_use,    total_mem_use,    total_mem_use2,    mem_use2 ));
end;

end
