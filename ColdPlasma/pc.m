function y=pc(x)
%------------------------------------------------------------
% Code for preconditioner
% compute M\r
% with M as the preconditioner matrix
%------------------------------------------------------------
global invM
% size(x)
% size(invM)
y = invM*x;




return;
