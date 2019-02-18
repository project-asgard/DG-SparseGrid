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
clear
close all

Lev = 2;
Deg = 2;
num_GridPoints = Deg * 2^Lev;
num_plot = Deg;


LInt = 0;
LEnd = 1;
Lmax = LEnd-LInt;
dx = Lmax/2^Lev;

CFL = .01;
dt = CFL*(dx)^2;
MaxT = ceil(1e-1/dt);

% Test 
ExactF = @(x,t)( exp(x) );
Source = @(x,t)( x-x );
BCFunc = @(x,t)( exp(x) );

ExactF = @(x,t)( sin(pi*x) );
Source = @(x,t)( x-x );
BCFunc = @(x,t)( sin(pi*x) );

Gam = 1;
f_bcL = 0; f_bcR = 0;
q_bcL = 1; q_bcR = 1;

Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,f_bcL,f_bcR); % equation for q
Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,@(x)1,@(x)0,q_bcL,q_bcR); % equation for f
Delta = Mat2*Mat1;
% MatBC = MatrixGradBC(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,q_bcL,q_bcR);

% Mat1 =   MatrixGradU(Lev,Deg,LInt,LEnd,-1,@(x)1,@(x)0,f_bcL,f_bcR); % equation for q
% Mat2 =   MatrixGradQ(Lev,Deg,LInt,LEnd,1,@(x)1,@(x)0,q_bcL,q_bcR); % equation for f
% Delta = Mat2*Mat1;

DoFs = (2^Lev*Deg);

II = speye(DoFs,DoFs);

Mat = kron(Delta,speye(DoFs,DoFs))+kron(speye(DoFs,DoFs),Delta);
% NewMat = speye(DoFs^2,DoFs^2) - dt*Mat;

F0 = zeros(DoFs,1);
time = 0;
n0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,time);
F0 = kron(n0,n0)*exp(-2*pi^2*time);

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);
[x_2D_plot,y_2D_plot] = meshgrid(x_node);
MM = kron(Meval,Meval);

bc = ComputeBC(Lev,Deg,LInt,LEnd,BCFunc,time,f_bcL,f_bcR);
bc = kron(bc,ones(DoFs,1)) +  kron(ones(DoFs,1),bc) ;


for T = 1 : 10%MaxT
    time = time + dt;
    bc0 = bc * exp(-2*pi^2*time);
    
    F1 = F0 + dt*(Mat)*F0 + dt*( bc0);
    val = exp(-2*pi^2*time)*sin(pi*x_2D_plot).*sin(pi*y_2D_plot);
    
%     F1 = NewMat\F0;
    mesh(x_2D_plot,y_2D_plot,reshape(MM*(F1),num_GridPoints,num_GridPoints));
%     hold on;
    F0 = F1;
    [norm(MM*(F1)-val(:)) max(abs(MM*(F1)-val(:)))]
%     val = exp(-2*pi^2*time)*sin(pi*x_2D_plot).*sin(pi*y_2D_plot);
    
%     mesh(x_2D_plot,y_2D_plot,val);
%     axis([0 1 0 1 0 1])
    title(['time = ', num2str(time)])
    hold off
    pause(0.1)
end

