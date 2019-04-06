function M = Helmholtz_1D(xMin,xMax,Lev,Deg,kappa,BC_opt)

sigma = 2;%1e-1i;%2;%1e-1i;%2;
rho = 0;%1e-1i;%0;%1e-3;%i;%0;%1e-1i;%4i;
Ord = Deg-1; % Max is Deg-1

% BCL = 'd';
% BCR = 'd';

if isempty(BC_opt.BCL)%~exist(BC_opt.BCL,'var') || isempty(BC_opt.BCL)
    BCL = 'N';
else
    BCL = BC_opt.BCL;
    FunBCL = BC_opt.FunBCL;
    FunExt = BC_opt.FunExt;
end
if isempty(BC_opt.BCR)%~exist('BC_opt.BCR','var') || isempty(BC_opt.BCR)
    BCR = 'N';
else
    BCR = BC_opt.BCR;
    FunBCR = BC_opt.FunBCR;
    FunExt = BC_opt.FunExt;    
end


FunRHS = @(x)(x-x-1);
FunExt = @(x)((1-cos(kappa*x)-sin(kappa)*sin(kappa*x)+1i*(1-cos(kappa))*sin(kappa*x))/kappa^2);
FunBCL = @(x)(x-x);
FunBCR = @(x)(x-x);

LMax = xMax - xMin;

N = 2^Lev; % # Grid Points
h = LMax/N;% size of the mesh
Jacobi = h/2;

% [Lev Deg kappa]% kappa^3*h^2/(Deg-1)^2]
%% Compute the 1D matrices
% compute the trace values for Polynomial and Derivatives
pL = legendre2(-1,Deg) * 1/sqrt(h);
pR = legendre2( 1,Deg) * 1/sqrt(h);
DDL = dlegendre2(-1,Deg,Deg-1) * 1/sqrt(h) * 2/h;
DDR = dlegendre2( 1,Deg,Deg-1) * 1/sqrt(h) * 2/h;
DDL = reshape(DDL,Deg,Deg)';
DDR = reshape(DDR,Deg,Deg)';

%---------------------------
% Matrices
%---------------------------
quad_num = 13; % 10 Gaussian Quadratures

[quad_x,quad_w]=lgwt(quad_num,-1,1);
pVal  = legendre2(quad_x,Deg)  * 1/sqrt(h);
tmp = dlegendre2(quad_x,Deg,1);
DpVal = tmp(:,:,2) * 1/sqrt(h) * 2/h;


DoF_1D = Deg*N;
A = sparse(DoF_1D,DoF_1D);
b = sparse(DoF_1D,1);
bc = sparse(DoF_1D,1);
uu = sparse(DoF_1D,1);

% Volume Integral
VolInt = [DpVal'*(quad_w.*DpVal)]*Jacobi;

% Trace Integral
% B:= <{u'},[v]> + <[u],{v'}>
% [B1 B2 B3]
B1 = (-pL)' * (DDR(2,:)/2) + (DDL(2,:)'/2) * ( pR);
B2L =( pR)' * (DDR(2,:)/2) + (DDR(2,:)'/2) * ( pR);
B2 = (-pL)' * (DDL(2,:)/2) + (DDL(2,:)'/2) * (-pL)+...
    ( pR)' * (DDR(2,:)/2) + (DDR(2,:)'/2) * ( pR);
B2R = (-pL)' * (DDL(2,:)/2) + (DDL(2,:)'/2) * (-pL);
B3 = ( pR)' * (DDL(2,:)/2) + (DDR(2,:)'/2) * (-pL);

% C:= I * p/h*<[u],[v]>
% [C1 C2 C3]
C1 = sigma*( (-pL)'*( pR) )* (Deg-1)/h;
C2L = sigma*( (pR)'*(pR) )* (Deg-1)/h;
C2 = sigma*( (-pL)'*(-pL)  + (pR)'*(pR) )* (Deg-1)/h;
C2R = sigma*( (-pL)'*(-pL) )* (Deg-1)/h;
C3 = sigma*( ( pR)'*(-pL) )* (Deg-1)/h;
% D:= I * sum_j (h/k)^(2*j-1)*<[D^ju],[D^jv]>
D1 = zeros(Deg,Deg);
D2 = zeros(Deg,Deg);
D3 = zeros(Deg,Deg);
D2L = zeros(Deg,Deg);
D2R = zeros(Deg,Deg);

for j = 1 : Ord%Deg-1
    coef = (h/(Deg-1))^(2*j-1);%*j;
    D1 = D1 + ( (-DDL(j+1,:))'*( DDR(j+1,:)) )*coef;
    D2 = D2 + ( (-DDL(j+1,:))'*(-DDL(j+1,:)) + DDR(j+1,:)'*DDR(j+1,:) )*coef;
    D2L = D2L + (   DDR(j+1,:)'*DDR(j+1,:) )*coef;
    D2R = D2R + ( (-DDL(j+1,:))'*(-DDL(j+1,:)) )*coef;
    D3 = D3 + ( ( DDR(j+1,:))'*(-DDL(j+1,:)) )*coef;
end


for Num = 0 : N-1
    xL = xMin + h*Num;
    xR = xL + h;
    
    % volume terms
    %     VolInt = VolInt;
    c = Deg*Num + [1:Deg];
    IndV = meshgrid(c);
    A = A + sparse(IndV',IndV,VolInt,DoF_1D,DoF_1D);
    
    if Num > 0 && Num < N-1
        % trace terms
        TraInt = [-B1+sigma*C1+rho*D1,...
                  -B2+sigma*C2+rho*D2,...
                  -B3+sigma*C3+rho*D3];
        
        IndU = [IndV-Deg,IndV,IndV+Deg];
        IndV = [IndV',IndV',IndV'];

    elseif Num == 0
        IndU = [IndV, IndV+Deg];
        IndV = [IndV', IndV'];
        if BCL == 'r' % robin bc
            TraInt = [-B2L+sigma*C2L+rho*D2L + 1i*kappa*(-pL)'*(-pL), ...
                      -B3+sigma*C3+rho*D3];
            
            b(c) = b(c) + pL'*FunBCL(xMin);
        elseif BCL == 'd' %
            TraInt = [-B2L+sigma*C2L+rho*D2L-(-pL)' * (DDL(2,:)/2)*2-(DDL(2,:)'/2) * (-pL)*2, ...
                      -B3+sigma*C3+rho*D3];
            
            b(c) = b(c) + (DDL(2,:)')*FunExt(xMin);
        end

    elseif Num == N-1        
        IndU = [IndV-Deg, IndV];
        IndV = [IndV', IndV'];
        if BCR == 'r' % robin bc
            TraInt = [-B1+sigma*C1+rho*D1, ...
                      -B2R+sigma*C2R+rho*D2R + 1i*kappa*(pR)'*(pR)];
            
            b(c) = b(c) + pR'*FunBCR(xMax);
        elseif BCR == 'd' %
            TraInt = [-B1+sigma*C1+rho*D1, ...
            -B2R+sigma*C2R+rho*D2R-(pR)'*(DDR(2,:)/2)*2- (DDR(2,:)'/2) * ( pR)*2];
            
            b(c) = b(c) - (DDR(2,:)')*FunExt(xMax);
        end
    end
        A = A + sparse(IndU,IndV,TraInt,DoF_1D,DoF_1D);

    % assemble the Rhs
    xi = xMin + h *(quad_x/2+1/2+Num);
    Val = pVal'*(quad_w.*FunRHS(xi))*Jacobi;
    b(c) = b(c)+Val;
    uu(c) = pVal'*(quad_w.*FunExt(xi))*Jacobi;
end

% FMWT_COMP = OperatorTwoScale(Deg,2^Lev);
% condest(A- kappa^2*speye(DoF_1D,DoF_1D))
% convert A to multiwavelet basis
% M = FMWT_COMP*A*FMWT_COMP';

% condest(B- kappa^2*speye(DoF_1D,DoF_1D))
% M = B - kappa^2*speye(DoF_1D,DoF_1D);
M = A;
end