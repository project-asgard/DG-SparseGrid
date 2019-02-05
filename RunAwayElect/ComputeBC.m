function bc = ComputeBC(Lev,Deg,LInt,LEnd,Fun,time,bcL,bcR)
%function ComputeBC to compute the bc term

%-----------------------------------------------------
L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

bc = sparse(DoF,1);

quad_num = 10;

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
% [quad_x,quad_w] = lgwt(quad_num,-1,1);
% p_val  = legendre(quad_x,Deg)  * 1/sqrt(h);
% 
% Jacobi = h/2;


WorkCel = 0;
if bcL == 0 % left cell is Dirichlet boundary
    c = [1:Deg];
    
    p_val  = legendre(-1,Deg)  * 1/sqrt(h);
    IntVal =  p_val'*(Fun(LInt,time)) ;
    
    bc(c) = - IntVal;
       
end

WorkCel = Tol_Cel_Num - 1;
if bcR == 0 % right cell is Dirichlet boundary
    
    c = Deg*WorkCel+[1:Deg];
    p_val  = legendre(1,Deg)  * 1/sqrt(h);
    IntVal =  p_val'*(Fun(LEnd,time));
    bc(c) = IntVal;
end


end