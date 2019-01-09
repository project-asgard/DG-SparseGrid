function Mat = MatrixMass(Lev,Deg,LInt,LEnd,FunCoef)
%function Mat to compute the Mass operator
% [FunCoef*f]
% Volum = - (FunCoef*f,v) 
%-----------------------------------------------------
if ~exist('FunCoef','var') || isempty(FunCoef)
    FunCoef = @(x)1;
end


L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

Mat = sparse(DoF,DoF);

quad_num = 10;



%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
[quad_x,quad_w] = lgwt(quad_num,-1,1);
p_val  = legendre(quad_x,Deg)  * 1/sqrt(h);

Jacobi = h/2;

for WorkCel = 0 : Tol_Cel_Num - 1
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    c = Deg*WorkCel+[1:Deg];
    
    xL = LInt + WorkCel*h;
    xR = xL + h;
    PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
    
    IntVal = [p_val'*(quad_w.*FunCoef(PhyQuad).*p_val)] * Jacobi;
    
    Mat = Mat + sparse(c'*ones(1,Deg),ones(Deg,1)*c,IntVal,DoF,DoF);
    
end

% figure;spy(Mat)

end