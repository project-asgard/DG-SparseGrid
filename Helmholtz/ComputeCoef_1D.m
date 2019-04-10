function [b] = ComputeCoef_1D(xMin,xMax,Lev,Deg,Fun)
%===============================================================
% This code computes the RHS terms for 1D source function f(x)
%===============================================================

DoF = 2^Lev * Deg;
b = sparse(DoF,1);


LMax = xMax - xMin;
N = 2^Lev; % # Grid Points
h = LMax/N;% size of the mesh
Jacobi = h/2;

quad_num = 10; % 10 Gaussian Quadratures
[quad_x,quad_w]=lgwt(quad_num,-1,1);
pVal  = legendre2(quad_x,Deg)  * 1/sqrt(h);

for Num = 0 : N-1
    xL = xMin + h*Num;
    xR = xL + h;
    
    Ind = Deg*Num + [1:Deg];
 
    % assemble the Rhs
    xi = xMin + h *(quad_x/2+1/2+Num);
    Val = pVal'*(quad_w.*Fun(xi))*Jacobi;
    
    b(Ind) = b(Ind)+Val;
end

end