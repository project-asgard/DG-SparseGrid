function [total_mem_use, mem_use] = mem_kron2( nr1, nc1, nr2, nc2, nvec )
% 
% [total_mem_use, mem_use] = mem_kron2( nr1, nc1, nr2, nc2, nvec )
% 
% estimate total amount of temporary space needed
%
idebug = 1;

%  Y =  (A2 * X) * transpose(A1)
%  computed as
%  step 1: Ytmp = A2 * X,   X is  nc2 by nc1 by nvec
%  step 2: Y = Ytmp * transpose(A1)
%
%  note Ytmp is  nr2 by nc1 by nvec
%
mem_use = nr2 * nc1 * nvec;

% ------------------------------------------
% no need for extra temporary array in kron1
% since it is just calling GEMM
% ------------------------------------------
total_mem_use = mem_use;

if (idebug >= 1),
  disp(sprintf('mem_kron2: mem_use=%g, total_mem_use=%g',
                           mem_use,    total_mem_use ));
end;

end
