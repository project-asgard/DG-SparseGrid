function bc = ComputeBC(Lev,Deg,xMin,xMax,FunL,FunR,time,bcL,bcR)
% function ComputeBC to compute the bc term
% This is the evaluation for two points on 1D
% Func*v|_xMin and Func*v|_xMax
%-----------------------------------------------------

L = xMax-xMin;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

bc = sparse(DoF,1);

p_L = legendre(-1,Deg) * 1/sqrt(h);
p_R = legendre(+1,Deg) * 1/sqrt(h);

WorkCel = 0;
if bcL == 'D' % if Dirichlet boundary?
    c = [1:Deg];
    IntVal =  p_L'*(FunL(xMin,time)) ; 
    bc(c) = - IntVal;   
end

WorkCel = Tol_Cel_Num - 1;
if bcR == 'D'    
    c = Deg*WorkCel+[1:Deg];
    IntVal =  p_R'*(FunR(xMax,time));
    bc(c) = IntVal;
end

end
