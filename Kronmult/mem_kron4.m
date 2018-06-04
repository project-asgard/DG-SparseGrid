
function [total_mem_use, mem_use] = mem_kron4( nr1, nc1, ...
                                               nr2, nc2, ...
                                               nr3, nc3, ...
                                               nr4, nc4, ...
                                               nvec )
%
% [total_mem_use, mem_use] = mem_kron4( nr1, nc1, nr2,nc2, nr3,nc3, nr4,nc4, nvec )
% 
% 
% estimate total amount of temporary space needed
%
idebug = 1;
%  Y =  { kron(A2,A3,A4) * X } * transpose(A1)
%
%  X is (nc4 by nc3 by nc2 by nc1) * nvec
%  computed as
%  step 1: Ytmp = kron3( A2,A3, A4, X )
%  step 2: Y = Ytmp * transpose(A1)
%
%  note Ytmp is  (nr4 * nr3 * nr2) * nvec3
%
nvec3 = nc1 * nvec;
mem_use = (nr4 * nr3 * nr2 ) * nvec3;

[total_mem_use3, mem_use3] = mem_kron3( nr2,nc2, nr3,nc3, nr4,nc4, nvec3 );

% ------------------------------------------
% no need for extra temporary array in kron1
% since it is just calling GEMM
% ------------------------------------------
total_mem_use = total_mem_use3 + mem_use;

if (idebug >= 1),
  disp(sprintf('mem_kron4: mem_use=%g, total_mem_use=%g, total_mem_use3=%g, mem_use3=%g', ...
                           mem_use,    total_mem_use,    total_mem_use3,    mem_use3 ));
end;



end
