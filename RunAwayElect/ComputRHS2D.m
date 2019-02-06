function rhs = ComputRHS2D(Lev,Deg,LInt,LEnd,Fun,time)
%function rhs to compute the rhs term
% This is the code to compute 2D problem
%-----------------------------------------------------
if ~exist('time','var') || isempty(time)
    time = 0;
end

L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF1D = Deg * Tol_Cel_Num;
DoF2D = DoF1D^2;

rhs = sparse(DoF2D,1);

quad_num = 10;


% compute the trace values
p_L = legendre(-1,Deg) * 1/sqrt(h);
p_R = legendre(+1,Deg) * 1/sqrt(h);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
[quad_x,quad_w] = lgwt(quad_num,-1,1);
p_val  = legendre(quad_x,Deg)  * 1/sqrt(h);
p_val2D = kron(p_val,p_val);
quad_w2D = kron(quad_w,quad_w);

Jacobi = h/2;
Jacobi2D = Jacobi^2;

for WorkCel_X = 0 : Tol_Cel_Num - 1
    c_x = Deg*WorkCel_X+[1:Deg];
    xL = LInt + WorkCel_X*h;
    xR = xL + h;
    PhyQuad_X = quad_x*(xR-xL)/2+(xR+xL)/2;
    PhyQuad_X = repmat(PhyQuad_X,1,10);
    PhyQuad_X = PhyQuad_X';
    PhyQuad_X = PhyQuad_X(:);
    
    for WorkCel_Y = 0 : Tol_Cel_Num-1
        c_y = Deg*WorkCel_Y+[1:Deg];
        yL = LInt + WorkCel_Y*h;
        yR = yL + h;
        PhyQuad_Y = quad_x*(yR-yL)/2+(yR+yL)/2;
        PhyQuad_Y = repmat(PhyQuad_Y,10,1);
        
        IntVal2D = p_val2D'*(quad_w2D.*Fun(PhyQuad_X,PhyQuad_Y,time)) * Jacobi2D;
        
%     c = Deg*WorkCel+[1:Deg];
%     
% 
%     
%     IntVal =  p_val'*(quad_w.*Fun(PhyQuad,time)) * Jacobi;
%         c = WorkCel_X * Tol_Cel_Num * Deg^2 + c_y;
        c = WorkCel_X * Deg * DoF1D + c_y'+[0:Deg-1]*DoF1D;
        c = c(:) ;
        rhs(c) = rhs(c) + IntVal2D;
    
    end
end

end