% Test 2D

% This is a test for 2D problem
% df/dt - Delta f = S
% Let q = grad f
% => df/dt - grad q = S
%    q = grad f

Lev = 5;
Deg = 2;
num_plot = 2;


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

ExactF = @(x,t)( 10*sin(pi*x) );
Source = @(x,t)( x-x );
BCFunc = @(x,t)( 10*sin(pi*x) );


% Source = sin(pi*x).*sin(pi*y)-2*pi*t*x.*cos(pi*x).*sin(pi*y)+2*pi^2*t*x.^2.*sin(pi*x).*sin(pi*y);

Gam = 1;
f_bcL = 0; f_bcR = 0;
q_bcL = 1; q_bcR = 1;

Mat1 = MatrixGrad(Lev,Deg,LInt,LEnd,-1,@(x)x.^2,@(x)0,f_bcL,f_bcR); % equation for q
Mat2 = MatrixGrad(Lev,Deg,LInt,LEnd, 1,@(x)1,@(x)0,q_bcL,q_bcR); % equation for f
Delta = Mat2*Mat1;


DoFs = (2^Lev*Deg);


Mat = kron(Delta,speye(DoFs,DoFs))+kron(speye(DoFs,DoFs),Delta);

F0 = zeros(DoFs,1);
time = 0;
n0 = ComputRHS(Lev,Deg,LInt,LEnd,ExactF,time);
F0 = kron(n0,n0)*exp(2*pi^2*time);

[x_node,Meval] = PlotDGData(Lev,Deg,LInt,LEnd,num_plot);
[x_2D_plot,y_2D_plot] = meshgrid(x_node);
MM = kron(Meval,Meval);

bc = ComputeBC(Lev,Deg,LInt,LEnd,BCFunc,time,f_bcL,f_bcR);
bc = kron(bc,bc);
for T = 1 : MaxT
    time = time + dt;
    bc0 = bc*exp(time);
    
    F1 = F0 + dt*(Mat)*F0 + dt*( bc0);
%     mesh(x_2D_plot,y_2D_plot,reshape(MM*(F0-F1),64,64),reshape(MM*(F0-F1),64,64))
    F0 = F1;
    
    val = 100*exp(2*pi^2*time)*sin(pi*x_2D_plot).*sin(pi*y_2D_plot);
    
    mesh(x_2D_plot,y_2D_plot,val);%reshape(MM*(F1),64,64))
    title(['time = ', num2str(time)])
    pause(0.1)
end

