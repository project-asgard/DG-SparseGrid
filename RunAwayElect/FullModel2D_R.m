% This is a test for FullModel 2D problem
% Gamma^C
% Term 1
% q = Ca * p^2 * df/dp
% Q = (Gp x I ) * F
% Term 2
% r = Cb * (1-x^2) * df/dx
% R = (I x Gx ) * F
% Term 3
% 1/p^2*d/dp [q] + 1/p^2 d/dp [Cf*p^2*f] + 1/p^4 d/dx[r]
% => 

Lev = 4;
Deg = 2;
num_plot = 2;

DoFs = (2^Lev*Deg);


LInt = -1;
LEnd = 1;

pInt = 0;
pEnd = 10;

Lmax = LEnd-LInt;
dx = Lmax/2^Lev;

CFL = 0.001;
dt = 1e-1;%2e-9;%CFL*(dx)^3;
MaxT = ceil(100/dt);


FullModel;

Gam = 1;
f_bcL = 1; f_bcR = 0;
q_bcL = 0; q_bcR = 1;

Mat_Mass_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(x.^2));

E = 0.0025;
Z = 1;
E = 0.0025;
tau = 10^5;
gamma = @(p)sqrt(1+(delta*p).^2);

% Term 1
% x*d/dp(p^2*E*f)

f_bcL = 1; f_bcR = 0;
q_bcL = 0; q_bcR = 1;

fp_bcL = 1; fp_bcR = 0;
qp_bcL = 0; qp_bcR = 1;

fx_bcL = 1; fx_bcR = 1;
qx_bcL = 0; qx_bcR = 0;


Mat_Term1_p = MatrixGrad(Lev,Deg,pInt,pEnd,1,@(x)x.^3.*gamma(x)/tau,@(x)0,fp_bcL,fp_bcR); 
Mat_Term1_x = MatrixMass(Lev,Deg,LInt,LEnd,@(x)(1-x.^2));
Mat_Term1 =  kron(Mat_Term1_p,Mat_Term1_x);

% Term 2
% r = 1/p^2 * d/dp*[p^2*Cf*f]
% R = (Gp * I ) * F
% Here we should use BC for x variable
Mat_Term2_p = MatrixMass(Lev,Deg,pInt,pEnd,@(x)(x.^2./(tau*gamma(x))));
Mat_Term2_x = MatrixGrad(Lev,Deg,pInt,pEnd,0,@(x)x.*(1-x.^2),@(x)0,fx_bcL,fx_bcR); 
Mat_Term2 = - kron(Mat_Term2_p,Mat_Term2_x);


Mat_Mass = kron(Mat_Mass_p,speye(DoFs,DoFs));
Inv = inv(Mat_Mass);

Mat_All = Mat_Term1 + Mat_Term2;


% Term 3
% 1/p^2*d/dp [q] + 1/p^2 d/dp [Cf*p^2*f] + 1/p^4 d/dx[r]
% => 

DoFs = (2^Lev*Deg)^2;

F0 = zeros(DoFs,1);

% Initial Condition
time = 0;
F0 = ComputRHS2D(Lev,Deg,[pInt LInt ],[pEnd LEnd ],Exa0,time);

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);
[p_node,Peval] = PlotDGData(Lev,Deg,pInt,pEnd,num_plot);
[x_2D_plot,y_2D_plot] = meshgrid(p_node,x_node);
MM = kron(Peval,Meval);
nz = size(x_2D_plot,1);

val_plot = reshape(MM*(F0),nz,nz);
surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
shading interp

Mat = Inv*Mat_All;
InvMat = inv(speye(DoFs,DoFs)-dt*Mat);

for T = 1 : MaxT
    
    rhs = F0-F0;%ComputRHS2D(Lev,Deg,LInt,LEnd,Source2D);
    
%     F1 = F0 + dt*(Mat)*F0 + dt*(Inv*rhs);
    
    % RK3
%     F1 = F0 + dt*(  Mat*F0 );
%     F2 = 3/4*F0+1/4*F1+1/4*dt*(Mat*F1);
%     Fn = 1/3*F0+2/3*F2+2/3*dt*(Mat*F2);
%     F0 = Fn;
    
    Fn = InvMat*F0;
    F0 = Fn;
%     mesh(x_2D_plot,y_2D_plot,reshape(MM*(F0-F1),64,64),reshape(MM*(F0-F1),64,64))
%     F0 = F1;
    
    time = time + dt;
%     bc0 = bc*exp(time);
    
%     val = exp(2*pi^2*time)*sin(pi*x_2D_plot).*sin(pi*y_2D_plot);
    
    val_plot = reshape(MM*(F0),nz,nz);
    surf(x_2D_plot,y_2D_plot,val_plot,val_plot)
    shading interp
%     view(0,90);
%     mesh(x_2D_plot,y_2D_plot,reshape(MM*(F1),64,64))
    title(['time = ', num2str(time)])
    pause(0.1)
end

