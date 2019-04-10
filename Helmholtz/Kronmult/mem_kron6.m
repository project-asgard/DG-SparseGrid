function [total_mem_use, mem_use] = mem_kron6( nr1, nc1, ...
                                               nr2, nc2, ...
                                               nr3, nc3, ...
                                               nr4, nc4, ...
                                               nr5, nc5, ...
                                               nr6, nc6, ...
                                               nvec )
%
% [total_mem_use, mem_use] = mem_kron6( nr1, nc1, nr2, nc2, ...
%                                       nr3, nc3, nr4, nc4, ...
%                                       nr5, nc5, nr6, nc6, nvec )
% 
% 
% estimate total amount of temporary space needed
%
idebug = 1;
%  Y =  { kron(A2,A3,A4,A5,A6) * X } * transpose(A1)
%
%  X is (nc6 by nc5 by nc4 by nc3 by nc2 by nc1) * nvec
%  computed as
%  step 1: Ytmp = kron5( A2,A3, A4, A5, A6, X )
%  step 2: Y = Ytmp * transpose(A1)
%
%  note Ytmp is  (nr6 * nr5 * nr4 * nr3 * nr2) * nvec5
%
nvec5 = nc1 * nvec;
mem_use = (nr6 * nr5 * nr4 * nr3 * nr2 ) * nvec5;

[total_mem_use5, mem_use5] = mem_kron5( nr2,nc2, nr3,nc3, ...
                                        nr4,nc4, nr5,nc5, ...
                                        nr6,nc6, nvec5 );

% ------------------------------------------
% no need for extra temporary array in kron1
% since it is just calling GEMM
% ------------------------------------------
total_mem_use = total_mem_use5 + mem_use;

if (idebug >= 1),
  disp(sprintf('mem_kron6: mem_use=%g, total_mem_use=%g, total_mem_use5=%g, mem_use5=%g', ...
                           mem_use,    total_mem_use,    total_mem_use5,    mem_use5 ));
end;

end
