function finit = Initial_Con_Grid(HASH,inverseHash,InCond)
%==========================================================
% This code is to generate initial conditions 
% which depend on Grids
% Input::
%       HASH and inverseHash
%       InCond
% Output::
%       finit
%==========================================================
dof = size(inverseHash,1);
finit = zeros(dof,1);

if Solver == 'VM'
    
elseif Solver == 'VP'
    
end