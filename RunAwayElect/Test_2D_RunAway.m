% Test 2D

% This is a test for 2D problem
% df/dt - Delta f = S
% Let q = grad f
% => df/dt - grad q = S
%    q = grad f

% Term 1
% q = Ca * p^2 * df/dp
% Q = (Gp x I ) * F
% Term 2
% r = Cb * (1-x^2) * df/dx
% R = (I x Gx ) * F
% Term 3
% 1/p^2*d/dp [q] + 1/p^2 d/dp [Cf*p^2*f] + 1/p^4 d/dx[r]
% => 

Lev = 3;
Deg = 2;
num_plot = 2;

DoFs = (2^Lev*Deg);


LInt = 0;
LEnd = 1;
Lmax = LEnd-LInt;
dx = Lmax/2^Lev;

CFL = 0.01;
dt = CFL*(dx)^3;
MaxT = ceil(1e-1/dt);

% Test 
ExactF = @(x,t)( exp(x) );
Source = @(x,t)( x-x );
BCFunc = @(x,t)( exp(x) );

ExactF = @(x,t)( sin(pi*x) );
Source = @(x,t)( x-x );
BCFunc = @(x,t)( sin(pi*x) );

ExactF = @(x,y,t)( sin(pi*x).*sin(pi*y)*exp(t) );


Ca = 1;
Cf = 0;
Cb = 0;

Source2D = @(x,y,t)x.^2.*( ...
    sin(pi*x).*sin(pi*y)...
    -Ca./x.^2.*exp(t).*( 2*pi*x.*cos(pi*x).*sin(pi*y)-pi^2*x.^2.*sin(pi*x).*sin(pi*y) )....
    -Cf./x.^2.*exp(t).*( 2*x.*sin(pi*x).*sin(pi*y)+pi*x.^2.*cos(pi*x).*sin(pi*y) )....
    -Cb./x.^4.*exp(t).*(-2*pi*y.*sin(pi*x).*cos(pi*y)-pi^2*(1-y.^2).*sin(pi*x).*sin(pi*y))...
    );


% Source = sin(pi*x).*sin(pi*y)-2*pi*t*x.*cos(pi*x).*sin(pi*y)+2*pi^2*t*x.^2.*sin(pi*x).*sin(pi*y);

Gam = 1;
f_bcL = 0; f_bcR = 0;
q_bcL = 1; q_bcR = 1;

Mat_Mass_p = MatrixMass(Lev,Deg,LInt,LEnd,@(x)(x.^2));

% Term 1
% q = Ca * p^2 * df/dp
% Q = (Gp x I ) * F

Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,@(x)x.^2,@(x)0,f_bcL,f_bcR); % equation for q
Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,q_bcL,q_bcR); % equation for f
Mat_Term1_p = Mat2*Mat1;

Mat_Term1 = kron(Mat_Term1_p,speye(DoFs,DoFs));

% Term 2
Mat_Term2_p = MatrixGrad(Lev,Deg,LInt,LEnd,1,@(x)x.^2,@(x)0,f_bcL,f_bcR); 

Mat_Term2 = kron(Mat_Term2_p,speye(DoFs,DoFs));

% Term 3
% r = Cb * (1-x^2) * df/dx
% R = (I x Gx ) * F
Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,@(x)(1-x.^2),@(x)0,f_bcL,f_bcR); % equation for q
Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,q_bcL,q_bcR); % equation for f
Mat_Term3_x = Mat2*Mat1;

Mat_Term3 = kron(Mat_Mass_p,Mat_Term3_x);

Mat_All = Ca*Mat_Term1 + Cf*Mat_Term2 + Cb*Mat_Term3;

Mat_Mass = kron(Mat_Mass_p,speye(DoFs,DoFs));
Inv = inv(Mat_Mass);

% All Terms
% rhs = ComputRHS2D(Lev,Deg,LInt,LEnd,Source2D);

% Term 3
% 1/p^2*d/dp [q] + 1/p^2 d/dp [Cf*p^2*f] + 1/p^4 d/dx[r]
% => 

DoFs = (2^Lev*Deg)^2;


% Mat = kron(Delta,speye(DoFs,DoFs))+kron(speye(DoFs,DoFs),Delta);

F0 = zeros(DoFs,1);
time = 0;
F0 = ComputRHS2D(Lev,Deg,LInt,LEnd,Exa0,time);
% F0 = kron(n0,n0)*exp(2*pi^2*time);

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);
[x_2D_plot,y_2D_plot] = meshgrid(x_node);
MM = kron(Meval,Meval);
nz = size(x_2D_plot,1);

val_plot = reshape(MM*(F0),nz,nz);
surf(x_2D_plot,y_2D_plot,val_plot,val_plot)

% val_plot = reshape(MM*(rhs),nz,nz);
% surf(x_2D_plot,y_2D_plot,val_plot,val_plot)

% bc = ComputeBC(Lev,Deg,LInt,LEnd,BCFunc,time,f_bcL,f_bcR);
% bc = kron(bc,bc);
for T = 1 : MaxT
    
    rhs = ComputRHS2D(Lev,Deg,LInt,LEnd,Source2D);
    
    F1 = F0 + dt*(Inv*Mat_All)*F0 + dt*(Inv*rhs);
%     mesh(x_2D_plot,y_2D_plot,reshape(MM*(F0-F1),64,64),reshape(MM*(F0-F1),64,64))
    F0 = F1;
    
    time = time + dt;
%     bc0 = bc*exp(time);
    
%     val = exp(2*pi^2*time)*sin(pi*x_2D_plot).*sin(pi*y_2D_plot);
    
    val_plot = reshape(MM*(F0),nz,nz);
    surf(x_2D_plot,y_2D_plot,val_plot,val_plot)

%     mesh(x_2D_plot,y_2D_plot,reshape(MM*(F1),64,64))
    title(['time = ', num2str(time)])
    pause(0.1)
end

