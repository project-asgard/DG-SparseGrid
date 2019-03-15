%% MATLAB (reference) version of the ASGarD solver
% New test for 2D Heat equation
% ToDo: The PDE definitions need more work
function  fk2dHeat
close all

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

pde = Diffusion2D;

gridType = 'SG';
nDims = 2;

for d=1:nDims
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde.dimensions{d}.deg,2^pde.dimensions{d}.lev);
end

deg = pde.dimensions{1}.deg;
lev = pde.dimensions{1}.lev;
Lmin = pde.dimensions{1}.domainMin;
Lmax = pde.dimensions{1}.domainMax;

BCFunc = @(x,t)(cos(pi*x)*exp(-2*pi^2*t));

% for time stepping
CFL = .01;
dx = 1/2^lev;
dt = CFL*(dx)^2;



% hash table
[HASH,HASHInv] = HashTable(lev,nDims,gridType);
nHash = numel(HASHInv);


t = 0;
% dealing with 1D operators
% df/dt = df^2/dx^2 + s
% assume s = 0
% let q = df/dx   (1)
% df/dt = dq/dx (2)
% (1) => Q = mat1*F + bc
% (2) => F = mat2*Q = mat2*(mat1*F + bc)
% here we use upwind for mat1 and downwind for mat2
% bc is Dirichlet boundary for variable f
% Denote Delta = mat2*mat1
Delta = coeff_matrix2(t,pde.dimensions{1},pde.term_1D);

dimension.lev = lev;
dimension.deg = deg;
dimension.domainMin = pde.dimensions{1}.domainMin;
dimension.domainMax = pde.dimensions{1}.domainMax;
dimension.FMWT = pde.dimensions{1}.FMWT;
dimension.BCL = 2; % Neumann
dimension.BCR = 2; % Neumann
term_1D.dat = [];
term_1D.LF = -1;       % Downwind Flux
term_1D.G = @(x,t,y)1; % Grad Operator
term_1D.type = 1;      % Grad Operator

% generate the matrix for variable q: This is only for handling boundary
% conditions
[mat2] = coeff_matrix2(t,dimension,term_1D);



% for initial condition
time = 0;
fval = initial_condition_vector(HASHInv,pde,time);

% boundary condition
% bc is the two points value for the first component
bc = ComputeBC(lev,deg,dimension.domainMin,dimension.domainMax,BCFunc,time,0,0);
bc = (pde.dimensions{1}.FMWT*bc);
% bc1 is the integration of \int_xMin^xMax f*v dx along the boundary 
% x = Const or y = Const
bc1 = ComputRHS(lev,deg,dimension.domainMin,dimension.domainMax,BCFunc,time);
bc1 = (pde.dimensions{1}.FMWT*bc1);

% construct the 2D boundary
ft = 1;
fList{1} = mat2*bc;
fList{2} = bc1;
bc3= combine_dimensions_D(fList,ft,HASHInv,pde);
fList{1} = bc1;
fList{2} = mat2*bc;
bc3= bc3+combine_dimensions_D(fList,ft,HASHInv,pde);


% 2D Matrix construction for sparse grids
DoFs = (2^lev*deg);
II = speye(DoFs,DoFs);

A_encode=GlobalMatrixSG_newCon(Delta,II,HASH,lev,deg,gridType);
B_encode=GlobalMatrixSG_newCon(II,Delta,HASH,lev,deg,gridType);
C_encode = [A_encode,B_encode];

[Meval_v,v_node,Meval_x,x_node]=matrix_plot(lev,lev,deg,...
    Lmin,Lmax, Lmin, Lmax,...
    pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);


for T = 1 : 100
    time = dt*T;
    
    % Euler time advance 
    % df^1 = df^0 + dt*A*f^0 + dt*s + dt*bc
    % just like Apply(A,f)
    ftmp = fval-fval; 
    for i=1:size(A_encode,2)
        
        tmpA=A_encode{i}.A1;
        tmpB=A_encode{i}.A2;
        IndexI=A_encode{i}.IndexI;
        IndexJ=A_encode{i}.IndexJ;
        
        
        nrA = size(tmpA,1);
        ncA = size(tmpA,2);
        nrB = size(tmpB,1);
        ncB = size(tmpB,2);
        
        ftmp(IndexI)=ftmp(IndexI) + ...
            reshape(tmpB * reshape(fval(IndexJ),ncB,ncA)*transpose(tmpA), nrB*nrA,1);
        
    end
    % bc term is cos(pi*x)*cos(pi*y)*exp(-2*pi^2*time)
    % thus we fully decomposite the x,y,t components
    ftmp = fval + dt*ftmp - dt* bc3 * exp(-2*pi^2*time);
    fval = ftmp;
    
    % plotting
    tmp = Multi_2D(Meval_v,Meval_x,fval,HASHInv,lev,deg);
    figure(1000)
    
    f2d0 = reshape(tmp,deg*2^lev,deg*2^lev)';
    
    [xx,vv]=meshgrid(x_node,v_node);
    
    ax1 = subplot(1,3,1);
    mesh(xx,vv,f2d0,'FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Lmin Lmax])
    view(-21,39)
    
    title(num2str(T))
    ax2 = subplot(1,3,2);
    val = cos(pi*xx).*cos(pi*vv)*exp(-2*pi^2*dt*T);
    mesh(xx,vv,val,'FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Lmin Lmax])
    view(-21,39)
    ax2 = subplot(1,3,3);
    mesh(xx,vv,val-f2d0,'FaceColor','interp','EdgeColor','none');
    axis([Lmin Lmax Lmin Lmax])
    title(num2str(max(abs(val(:)-f2d0(:)))))
    view(-21,39)
    pause(0.1)
end







end

