function bc = compute_boundary_condition(pde,g_func,dV,t,Lev,Deg,xMin,xMax,bc_func,LorR)
% function bc = ComputeBC(Lev,Deg,xMin,xMax,Fun,time,LorR)
% function ComputeBC to compute the bc term
% This is the evaluation for two points on 1D
% Func*v|_xMin and Func*v|_xMax
%-----------------------------------------------------

p = pde.params;

L = xMax-xMin;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
small_dx = h*1e-7; %Parameter to fix problem with Inhomogeneous BC
DoF = Deg * Tol_Cel_Num;

bc = sparse(DoF,1);

p_L = lin_legendre(-1,Deg) * 1/sqrt(h); % TODO : this happens in multiple places. Consolidate. 
p_R = lin_legendre(+1,Deg) * 1/sqrt(h);

if strcmp(LorR,'L')
    
    WorkCel = 0;
    c = [1:Deg];
    IntVal =  p_L'*bc_func(xMin,p,t);
    if isfinite(g_func(xMin,p,t)) %make sure g_func is finite
        
        bc(c) = -g_func(xMin,p,t).*IntVal.*dV(xMin,p,t);
    
    else
        
        xLclose = xMin + small_dx; %Move away from xMin by small_dx
        bc(c) = -g_func(xLclose,p,t).*IntVal.*dV(xLclose,p,t);
        
    end
else
    
    WorkCel = Tol_Cel_Num - 1;
    c = Deg*WorkCel+[1:Deg];
    IntVal =  p_R'*bc_func(xMax,p,t);
    
    if isfinite(g_func(xMax,p,t)) %make sure g_func is finite
        
        bc(c) = g_func(xMax,p,t).*IntVal.*dV(xMax,p,t);
        
    else
        
        xRclose = xMax - small_dx; %Move away from xMax by small_dx
        bc(c) = g_func(xRclose,p,t).*Intval.*dV(xRclose,p,t);
        
    end
    
end

end
