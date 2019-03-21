function rhsList = ComputRHS(Lev,Deg,LInt,LEnd,Fun,time)
%function rhs to compute the rhs term

%-----------------------------------------------------
if ~exist('time','var') || isempty(time)
    time = 0;
end

L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

quad_num = 10;

% compute the trace values
% p_L = legendre(-1,Deg) * 1/sqrt(h);
% p_R = legendre(+1,Deg) * 1/sqrt(h);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
[quad_x,quad_w] = lgwt(quad_num,-1,1);
p_val  = legendre(quad_x,Deg)  * 1/sqrt(h);

Jacobi = h/2;


% rhs = sparse(DoF,1);
% for WorkCel = 0 : Tol_Cel_Num - 1
%     %---------------------------------------------
%     % (funcCoef*q,d/dx p)
%     %---------------------------------------------
%     c = Deg*WorkCel+[1:Deg];
%     
%     xL = LInt + WorkCel*h;
%     xR = xL + h;
%     PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
%     
%     IntVal =  p_val'*(quad_w.*Fun{1}(PhyQuad,time)) * Jacobi;
% 
%     rhs(c) = rhs(c) + IntVal;
%     
% end

nDims = numel(Fun);
for d = 1 : nDims
    rhs = sparse(DoF,1);
    for WorkCel = 0 : Tol_Cel_Num - 1
        %---------------------------------------------
        % (funcCoef*q,d/dx p)
        %---------------------------------------------
        c = Deg*WorkCel+[1:Deg];
        
        xL = LInt + WorkCel*h;
        xR = xL + h;
        PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
        
        IntVal =  p_val'*(quad_w.*Fun{d}(PhyQuad,time)) * Jacobi;
        
        rhs(c) = rhs(c) + IntVal;
        
    end
    rhsList{d} = rhs;
end

end
