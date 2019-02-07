function Mat = MatrixBurgers(Lev,Deg,LInt,LEnd,FluxVal,FunCoef,bcL,bcR)
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

Mat = sparse(DoF,DoF);

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

for WorkCel = 0 : Tol_Cel_Num - 1
    %---------------------------------------------
    % (funcCoef*q,d/dx p)
    %---------------------------------------------
    c = Deg*WorkCel+[1:Deg];
    
    xL = LInt + WorkCel*h;
    xR = xL + h;
    PhyQuad = quad_x*(xR-xL)/2+(xR+xL)/2;
    
    ProjF = p_val*FunCoef(c); % solution on each quadratures
    FunVal = p_val'*(quad_w.*ProjF.^2) * Jacobi;
    
%     FunVal = legendre( 0,Deg)* 1/sqrt(h)*FunCoef(c);
    if WorkCel > 0
        ProjF = p_val*FunCoef(c-Deg); % solution on each quadratures
        FunValPre = p_val'*(quad_w.*ProjF.^2) * Jacobi;
%         FunValPre = legendre( 0,Deg)* 1/sqrt(h)*FunCoef(c-Deg);       
    else
        FunValPre = FunVal;
    end
    if WorkCel < Tol_Cel_Num - 1
        ProjF = p_val*FunCoef(c+Deg); % solution on each quadratures
        FunValLat = p_val'*(quad_w.*ProjF.^2) * Jacobi;
%         FunValLat = legendre( 0,Deg)* 1/sqrt(h)*FunCoef(c+Deg);       
    else
        FunValLat = FunVal;
    end  
    
    IntVal = ...
        - [Dp_val'*(quad_w.*FunVal.*p_val)] * Jacobi;
    
    Mat = Mat + sparse(c'*ones(1,Deg),ones(Deg,1)*c,IntVal,DoF,DoF);
    %----------------------------------------------
    % -<funcCoef*{q},p>
    %----------------------------------------------
    % Numerical Flux is defined as
    % Flux = {{f}} + C/2*[[u]]
    %      = ( f_L + f_R )/2 + FunCoef*( u_R - u_L )/2
    % [[v]] = v_R - v_L

%     TraVal = [...
%         (-p_L)' * FunCoef(xL) * p_R/2 + FluxVal * abs(FunCoef(xL))/2 * (-p_L)' *   p_R,...
%         (-p_L)' * FunCoef(xL) * p_L/2 + FluxVal * abs(FunCoef(xL))/2 * (-p_L)' * (-p_L),...% xL
%         ( p_R)' * FunCoef(xR) * p_R/2 + FluxVal * abs(FunCoef(xR))/2 * ( p_R)' *   p_R,...
%         ( p_R)' * FunCoef(xR) * p_L/2 + FluxVal * abs(FunCoef(xR))/2 * ( p_R)' * (-p_L),...% xR
%         ];
%     if bcL == 0 && WorkCel == 0
%         TraVal = [...
%             (-p_L)' * (p_L-p_L),...
%             (-p_L)' * (p_L-p_L),...% xL
%             ( p_R)' * FunCoef(xR) * p_R/2 + FluxVal * abs(FunCoef(xR))/2 * ( p_R)' *   p_R,...
%             ( p_R)' * FunCoef(xR) * p_L/2 + FluxVal * abs(FunCoef(xR))/2 * ( p_R)' * (-p_L),...% xR
%             ];
% 
%     end
%     if bcR == 0 && WorkCel == Tol_Cel_Num - 1
%         TraVal = [...
%             (-p_L)' * FunCoef(xL) * p_R/2 + FluxVal * abs(FunCoef(xL))/2 * (-p_L)' *   p_R,...
%             (-p_L)' * FunCoef(xL) * p_L/2 + FluxVal * abs(FunCoef(xL))/2 * (-p_L)' * (-p_L),...% xL
%             ( p_R)' * (p_R-p_R),...
%             ( p_R)' * (p_R-p_R),...% xR
%             ];
%     end
    
	TraVal = [...
        (-p_L)' * FunValPre * p_R/2 + FluxVal * abs(FunValPre)/2 * (-p_L)' *   p_R,...
        (-p_L)' * FunVal * p_L/2 + FluxVal * abs(FunVal)/2 * (-p_L)' * (-p_L),...% xL
        ( p_R)' * FunVal * p_R/2 + FluxVal * abs(FunVal)/2 * ( p_R)' *   p_R,...
        ( p_R)' * FunValLat * p_L/2 + FluxVal * abs(FunValLat)/2 * ( p_R)' * (-p_L),...% xR
        ];
    if bcL == 0 && WorkCel == 0
        TraVal = [...
            (-p_L)' * (p_L-p_L),...
            (-p_L)' * (p_L-p_L),...% xL
            ( p_R)' * FunVal * p_R/2 + FluxVal * abs(FunVal)/2 * ( p_R)' *   p_R,...
            ( p_R)' * FunValLat * p_L/2 + FluxVal * abs(FunValLat)/2 * ( p_R)' * (-p_L),...% xR
            ];

    end
    if bcR == 0 && WorkCel == Tol_Cel_Num - 1
        TraVal = [...
            (-p_L)' * FunValPre * p_R/2 + FluxVal * abs(FunValPre)/2 * (-p_L)' *   p_R,...
            (-p_L)' * FunVal * p_L/2 + FluxVal * abs(FunVal)/2 * (-p_L)' * (-p_L),...% xL
            ( p_R)' * (p_R-p_R),...
            ( p_R)' * (p_R-p_R),...% xR
            ];
    end
    % Adding trace value to matrix
    RowInd = [c'*ones(1,Deg) c'*ones(1,Deg) c'*ones(1,Deg) c'*ones(1,Deg)];
    ColInd = [ones(Deg,1)*(c-Deg),ones(Deg,1)*c,ones(Deg,1)*c,ones(Deg,1)*(c+Deg)];
    
    
    if WorkCel == 0
        Iu = RowInd(:,Deg+1:end);
        Iv = ColInd(:,Deg+1:end);
        Val = TraVal(:,Deg+1:end);
    elseif WorkCel == Tol_Cel_Num - 1
        Iu = RowInd(:,1:3*Deg);
        Iv = ColInd(:,1:3*Deg);
        Val = TraVal(:,1:3*Deg);
    else
        Iu = RowInd;
        Iv = ColInd;
        Val = TraVal;
    end
    
    
    
    Mat = Mat + sparse(Iu,Iv,Val,DoF,DoF);
    
  
    
    
end

% figure;spy(Mat)

end