function [total_mem_use, mem_use] = mem_kron5( nr1, nc1, ...
                                               nr2, nc2, ...
                                               nr3, nc3, ...
                                               nr4, nc4, ...
                                               nr5, nc5, ...
                                               nvec )
%
% [total_mem_use, mem_use] = mem_kron5( nr1, nc1, nr2, nc2, ...
%                                       nr3, nc3, nr4, nc4, ...
%                                       nr5, nc5, nvec )
% 
% 
% estimate total amount of temporary space needed
%
idebug = 1;
%  Y =  { kron(A2,A3,A4,A5) * X } * transpose(A1)
%
%  X is (nc5 by nc4 by nc3 by nc2 by nc1) * nvec
%  computed as
%  step 1: Ytmp = kron3( A2,A3, A4, X )
%  step 2: Y = Ytmp * transpose(A1)
%
%  note Ytmp is  (nr5 * nr4 * nr3 * nr2) * nvec4
%
nvec4 = nc1 * nvec;
mem_use = (nr5 * nr4 * nr3 * nr2 ) * nvec4;

[total_mem_use4, mem_use4] = mem_kron4( nr2,nc2, nr3,nc3, ...
                                        nr4,nc4, nr5,nc5, nvec4 );

% ------------------------------------------
% no need for extra temporary array in kron1
% since it is just calling GEMM
% ------------------------------------------
total_mem_use = total_mem_use4 + mem_use;

if (idebug >= 1),
  disp(sprintf('mem_kron5: mem_use=%g, total_mem_use=%g, total_mem_use4=%g, mem_use4=%g', ...
                           mem_use,    total_mem_use,    total_mem_use4,    mem_use4 ));
end;

end
