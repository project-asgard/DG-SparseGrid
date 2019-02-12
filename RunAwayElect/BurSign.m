function [uL,uR] = BurSign(LInt,LEnd,Lev,Deg,uold)
L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;


% compute the trace values
p_L = legendre(-1,Deg) * 1/sqrt(h);
p_R = legendre(+1,Deg) * 1/sqrt(h);

uL = p_L*uold(1:Deg);
uR = p_R*uold(end-Deg+1:end);