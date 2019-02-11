function [FF,ProjU] = VectorBurgers(Lev,Deg,LInt,LEnd,FluxVal,uold,bcL,bcR)
%function Mat to compute the Grad operator
% The FunCoef is a nonlinear term in the PDEs
% FunCoef*d/dx[f]
% Trace = FunCoef*f|_R - FunCoef*f|_L+
% Volum = - (FunCoef*f,d/dx v)
% Flux = {f}+FluxVal*abs(FunCoef)*[f]/2
% FluxVal ::
% FluxVal = 0 --> Central Flux
% FluxVal = 1 --> Upwind Flux
%
% FunCoef
% FunCoef2
%
% bcL and bcR ::
% bc = 0 --> Dirichlet; bc = 1 --> Neumann
%-----------------------------------------------------
if ~exist('FluxVal','var') || isempty(FluxVal)
    FluxVal = 0;
end

if ~exist('FunCoef','var') || isempty(FunCoef)
    FunCoef = @(x)1;
end


% bcL = 0 --> Dirichlet; bcL = 1 --> Neumann
% bcR = 0 --> Dirichlet; bcR = 1 --> Neumann
% The default value is 1 as Neumann boundary condition
if ~exist('bcL','var') || isempty(bcL)
    bcL = 1;
end

if ~exist('bcR','var') || isempty(bcR)
    bcR = 1;
end

L = LEnd-LInt;
Tol_Cel_Num = 2^(Lev);
h = L  / Tol_Cel_Num;
DoF = Deg * Tol_Cel_Num;

FF = sparse(DoF,1);

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
Dp_val = dlegendre(quad_x,Deg) * 1/sqrt(h) * 2/h;

Jacobi = h/2;

% First compute the projection of u^2/2 on each cell
MaxC = 0;
for WorkCel = 0 : Tol_Cel_Num - 1
    c = Deg*WorkCel+[1:Deg];
    
    xL = LInt + WorkCel*h;
    xR = xL + h;
    PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
    
    EleFVal = uold(c); % solution on each quadratures
    NonVal = zeros(quad_num,1);
    
    MaxC = max([MaxC;abs(p_val*EleFVal)]);
    
    for i = 1 : Deg
        for j = 1 : Deg
            NonVal = NonVal + (p_val(:,i)*EleFVal(i)).*(p_val(:,j)*EleFVal(j));
        end
    end
    
    ProjU(c,1) = p_val'*(quad_w.*NonVal/2)* Jacobi;
    
end

MaxC

for WorkCel = 0 : Tol_Cel_Num - 1
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    c = Deg*WorkCel+[1:Deg];
    
    xL = LInt + WorkCel*h;
    xR = xL + h;
    PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
    
    EleFVal = p_val*ProjU(c); % solution on each quadratures


    
%     FunVal = legendre( 0,Deg)* 1/sqrt(h)*FunCoef(c);
    if WorkCel > 0
        TraceFL = (p_R*ProjU(c-Deg) + p_L*ProjU(c))/2 +... % avg = {u^2/2}
            + MaxC/2 * (p_R*ProjU(c-Deg) - p_L*ProjU(c) ); % Jum = [u^2/2] solution on each quadratures
    else
        TraceFL = ProjU(c)-ProjU(c);%(p_L*ProjU(c));
        %ProjU(c)-ProjU(c);%p_L*ProjU(c)/2 - MaxC/2 * p_L*ProjU(c);%- p_L*ProjU(c);
    end
    if WorkCel < Tol_Cel_Num - 1
        TraceFR = (p_R*ProjU(c) + p_L*ProjU(c+Deg))/2 + ...
            + MaxC/2 * (p_R*ProjU(c) - p_L*ProjU(c+Deg));       
    else
        TraceFR = (p_R*ProjU(c) );
        %ProjU(c)-ProjU(c);%p_R*ProjU(c)/2 + MaxC/2 * p_R*ProjU(c);%p_R*ProjU(c);
    end  
    
    IntVal = ...
        - [Dp_val'*(quad_w.*EleFVal)] * Jacobi;
    
    FF = FF + sparse(c,ones(Deg,1),IntVal,DoF,1);
    
    TraceVal = - p_L * TraceFL + p_R * TraceFR;
    FF = FF + sparse(c,ones(Deg,1),TraceVal,DoF,1);
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    % Numerical Flux is defined as
    % Flux = {{f}} + C/2*[[u]]
    %      = ( f_L + f_R )/2 + FunCoef*( u_R - u_L )/2
    % [[v]] = v_R - v_L

    
end

% figure;spy(Mat)

end