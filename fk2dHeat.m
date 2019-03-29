%% MATLAB (reference) version of the ASGarD solver
% New test for 2D Heat equation
% ToDo: The PDE definitions need more work
function  fk2dHeat
close all

format short e
folder = fileparts(which(mfilename));
addpath(genpath(folder));

pde = diffusion2;

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
% but the Boundary Function is defined inside coeff_matrix2

%% Boundary conditions
% bc is the two points value for the first component
% bc1 is the integration of \int_xMin^xMax f*v dx along the boundary
% x = Const or y = Const

time = 0;
for d = 1:nDims
    
    [Delta,bcL_tmp,bcR_tmp] = coeff_matrix2(t,pde.dimensions{d},pde.terms{1}{1},d);
    bcL{d} = bcL_tmp;
    bcR{d} = bcR_tmp;
    
    dim = pde.dimensions{d};
    
    bcL1{d} = ComputeRHS(nDims,dim,dim.BCL_fList);
    bcR1{d} = ComputeRHS(nDims,dim,dim.BCR_fList);
    
end


%%
% Construct the RHS piece of the Neumann BCs
%
% TODO
%

% for initial condition
fval = initial_condition_vector(HASHInv,pde,time);

% 2D Matrix construction for sparse grids
DoFs = (2^lev*deg);
II = speye(DoFs,DoFs);

A_encode=GlobalMatrixSG_newCon(Delta,II,HASH,lev,deg,gridType);
B_encode=GlobalMatrixSG_newCon(II,Delta,HASH,lev,deg,gridType);

A_encode = [A_encode,B_encode];

[Meval_v,v_node,Meval_x,x_node]=matrix_plot(lev,lev,deg,...
    Lmin,Lmax, Lmin, Lmax,...
    pde.dimensions{1}.FMWT,pde.dimensions{2}.FMWT);


for T = 1 : 100
    time = dt*T;
    
    %%
    % Construct the time-dependent RHS piece of the Dirichlet BCs
    
    bc3 = zeros(deg^nDims*nHash,1);
    for d=1:nDims
        
        dim = pde.dimensions{d};
        
        ftL = dim.BCL_fList{nDims+1}(time)
        ftR = dim.BCR_fList{nDims+1}(time);
        
        %%
        % For each dimensional boundary
        clear fListL fListR;
        for d2=1:nDims
            
            if d == d2
                fListL{d2} = bcL{d};
                fListR{d2} = bcR{d};
            else
                fListL{d2} = bcL1{d}{d2};
                fListR{d2} = bcR1{d}{d2};
            end
        end
        
        bc3 = bc3 + combine_dimensions_D(fListL,ftL,HASHInv,pde);
        bc3 = bc3 + combine_dimensions_D(fListR,ftR,HASHInv,pde);
        
    end
    
    % Euler time advance for zero source
    % df^1 = df^0 + dt*A*f^0  + dt*bc
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
    ftmp = fval + dt*ftmp - dt* bc3;% * exp(-2*pi^2*time);
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
