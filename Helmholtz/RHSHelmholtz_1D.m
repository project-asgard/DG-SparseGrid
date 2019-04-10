function [b,bc] = RHSHelmholtz_1D(xMin,xMax,Lev,Deg,FunRHS,BC_opt)
%===============================================================
% This code computes the RHS terms for 1D source function f(x)
%===============================================================

if isempty(BC_opt.BCL)%~exist('BC_opt.BCL','var') || isempty(BC_opt.BCL)
    BCL = 'N';
else
    BCL = BC_opt.BCL;
    FunBCL = BC_opt.FunBCL;
    FunExt = BC_opt.FunExt;
end
if isempty(BC_opt.BCR)%~exist('BC_opt.BCR','var') || isempty(BC_opt.BCR)
    BCR = 'N';
else
    BCR = BC_opt.BCR;
    FunBCR = BC_opt.FunBCR;
    FunExt = BC_opt.FunExt;    
end

DoF = 2^Lev * Deg;
b = sparse(DoF,1);
bc = sparse(DoF,1);

LMax = xMax - xMin;
N = 2^Lev; % # Grid Points
h = LMax/N;% size of the mesh
Jacobi = h/2;

quad_num = 10; % 10 Gaussian Quadratures
[quad_x,quad_w]=lgwt(quad_num,-1,1);
pVal  = legendre2(quad_x,Deg)  * 1/sqrt(h);

pL = legendre2(-1,Deg) * 1/sqrt(h);
pR = legendre2( 1,Deg) * 1/sqrt(h);

DDL = dlegendre2(-1,Deg,Deg-1) * 1/sqrt(h) * 2/h;
DDR = dlegendre2( 1,Deg,Deg-1) * 1/sqrt(h) * 2/h;
DDL = reshape(DDL,Deg,Deg)';
DDR = reshape(DDR,Deg,Deg)';

for Num = 0 : N-1
    xL = xMin + h*Num;
    xR = xL + h;
    
    Ind = Deg*Num + [1:Deg];
    
    if Num == 0
        
        if BCL == 'r' % Robin B.C.
            
            bc(Ind) = bc(Ind) + pL'*FunBCL(xMin);
            
        elseif BCL == 'd' % Dirichlet B.C.
            
            bc(Ind) = bc(Ind) + (DDL(2,:)')*FunExt(xMin);
            
        end
    end
    
    if Num == N-1
        
        if BCR == 'r' % Robin B.C.
            
            bc(Ind) = bc(Ind) + pR'*FunBCR(xMax);
            
        elseif BCR == 'd' % Dirichlet B.C.
            
            bc(Ind) = bc(Ind) - (DDR(2,:)')*FunExt(xMax);
            
        end
    end
    
    
    % assemble the Rhs
    xi = xMin + h *(quad_x/2+1/2+Num);
    Val = pVal'*(quad_w.*FunRHS(xi))*Jacobi;
    
    b(Ind) = b(Ind)+Val;
end

end