function bc = ComputeBC(Lev,Deg,LInt,LEnd,Fun,time,bcL,bcR,FluxVal,FunCoef,MaxC)
%function ComputeBC to compute the bc term

%-----------------------------------------------------
if ~exist('FluxVal','var') || isempty(FluxVal)
    FluxVal = 0;
end
if ~exist('FunCoef','var') || isempty(FunCoef)
    FunCoef = @(x)(x-x+1);
end
if ~exist('MaxC','var') || isempty(MaxC)
    MaxC = 1;
end

L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

bc = sparse(DoF,1);


%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
% [quad_x,quad_w] = lgwt(quad_num,-1,1);
% p_val  = legendre(quad_x,Deg)  * 1/sqrt(h);
% 
% Jacobi = h/2;

p_L = legendre(-1,Deg) * 1/sqrt(h);
p_R = legendre(+1,Deg) * 1/sqrt(h);

WorkCel = 0;

if bcL == 0 
    c = [1:Deg];
    p_val  = legendre(-1,Deg)  * 1/sqrt(h);
    IntVal =  p_val'*(Fun(LInt,time)) ; 
    bc(c) = - IntVal;   
end

WorkCel = Tol_Cel_Num - 1;

if bcR == 0    
    c = Deg*WorkCel+[1:Deg];
    p_val  = legendre(1,Deg)  * 1/sqrt(h);
    IntVal =  p_val'*(Fun(LEnd,time));
    bc(c) = IntVal;
end


end
