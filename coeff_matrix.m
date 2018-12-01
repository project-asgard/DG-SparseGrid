%% Construct 1D coefficient matrix
% This routine returns a 2D array representing an operator coefficient
% matrix for a single-D. Each term in a PDE requires D many coefficient
% matricies. These operators can only use the supported types below.
%
function coeff_matrix(t,lev,deg,type,xMin,xMax,BCL,BCR,LF,FMWT)

%% Inputs
% lev  : number of levels in hierachical basis
% deg  : order of legendre basis
% type : object of type term (see below)
% xMin : minimum of this D domain range
% xMax : maximum of this D domain range
% BCL  : Boundary Condition (Left),  select from list of BCs below
% BCR  : Boundary Condition (Right), select from list of BCs below
% LF   : =0 for CF flux, or =C for LF flux with global C value
% FMWT : Forward Multi-Wavelet Transform matrix

%% "type" object
% This object contains the information to specify the type of the operator.
% It defines the following ...
%  coeff_type : int from coeff_type list below
%  TD         : 0 or 1 (is this term time dependent)
%  g1(x,t)    : function handle to g1 function
%  g2(x,t)    : function handle to g2 function
%  g3(x,t)    : function handle to g3 function

%% Boundary condition types (BCL and BCR)
% 0 == periodic
% 1 == dirichlet (set value of solution)
% 2 == neumann   (set first derivative of solution)

%% Available coefficient types (coeff_type)
% 1 == FuncGrad
% 2 == FuncMass
% 3 == ???

%% Note on global vs local Lax-Friedrichs (LF) flux
% We do not (cannot) use local upwinding or LF because selecting
% either the sign of the flow field or the value of the coefficient C could
% be multivalued within the multi-D solution for a single-D coeff_matrix.

%function [vMassV,GradV,GradX,DeltaX,FluxX,FluxV]=matrix_coeff_TI(Lev_x,lev,k,Lmin,Lmax,Vmin,Vmax,FMWT_COMP_x,FMWT_COMP_v)
%=============================================================
% Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:  vMassV:  int_v v * l_i(v)  * l_j(v) dv
%               GradV :  int_v 1 * l_i(v)' * l_j(v) dv
% Poisson Solver:
%   Operators: DeltaX: int_x (m_i(x))''*m_j(x)dx
%   This equation is solved by LDG methods
% Maxwell Solver: (Missing)
%   Operators: CurlX: int_x curl(m_i(x))*m_j(x)dx
% Input: Lev, k, dim, Lmax, Vmax
% P.S. This is the full-grid version
%=============================================================

%% TODO ...
% * Choice of flux (may require input C)
% * Other BCs
% * Picking which term type

%%
% Set number of quatrature points (should this be order dependant?)
quad_num = 10;

%%
%  Get quadrature points and weights.
%  quad_x(:) is quad_num by 1
%  quad_w(:) is quad_num by 1
[quad_x,quad_w] = lgwt(quad_num,-1,1);

%%
%  Compute the trace values (values at the left and right of each element for all k)
%  p_1(:) is 1 by deg
%  p_2(:) is 1 by deg
p_1 = legendre(-1,deg);
p_2 = legendre(+1,deg);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
p_val  = legendre(quad_x,deg);
Dp_val = dlegendre(quad_x,deg);

%%
% Setup jacobi of variable x and define coeff_mat
N = 2^(lev);
h = (Vmax-Vmin) / N;
jacobi = h;
dof_1D = deg * N;

vMassV = sparse(dof_1D,dof_1D);
GradV  = sparse(dof_1D,dof_1D);
FluxV  = sparse(dof_1D,dof_1D);

%% Loop over all elements in this D
for i=0:N-1
    
    %%
    % Get index ranges for ...
    
    %%
    %  Current element
    c1 = deg*i+1;
    c2 = deg*(i+1);
    c = c1:c2;
    
    %%
    % Previous element
    p1 = deg*(i-1)+1;
    p2 = deg*i;
    p = p1:p2;
    
    %%
    % Later element
    l1 = deg*(i+1)+1;
    l2 = deg*(i+2);
    l = l1:l2;
    
    %%
    % Map from [-1,1] to physical domain
    
    x = (((quad_x+1)/2+i)*h+Vmin);
    
    %%
    % Build Average (AVG) and Jump (JMP) operators
    
    val_AVG=(1/h)*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
        p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    val_JMP=1/h*[p_1'*p_2 -p_1'*p_1,... % for x1
        -p_2'*p_2  p_2'*p_1]/2; % for x2 %%
    
    
    %%
    % Perform volume integral
    val_vMassV = p_val'*(p_val.*x.*quad_w)*jacobi/2/h;
    val_gradV  = 1/h*[Dp_val'*(quad_w.*p_val)];
    
    Iu=meshgrid(k*i+1:k*(i+1));
    
    vMassV=vMassV+sparse(Iu',Iu,val_vMassV,dof_1D,dof_1D);
    GradV =GradV +sparse(Iu',Iu,val_AVG,dof_1D,dof_1D);
    
    
    %%
    % Setup numerical flux choice (interior elements only)
    
    if i<N-1 && i>0
        
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
        Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        
    end
    
    %%
    % Setup boundary conditions
        
    if strcmp(BC,'periodic')
        
        %%
        % Left boundary
        if i==0
            Iu=[meshgrid([k*(N-1)+1:k*(N)]),meshgrid(c),meshgrid(c),meshgrid(l)];
            Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        end
        
        %%
        % Right boundary
        if i==N-1
            Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
            Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
        end
        
    end
    
    if strcmp(BC,'dirichlet')
        
        %%
        % TODO
        
    end
    
    if strcmp(BC,'neumann')
        
        %%
        % TODO
        
    end
    
    %%
    % Apply flux choice / BCs
    
    GradV = GradV - sparse(Iv,Iu,val_AVG,dof_1D,dof_1D);
    FluxV = FluxV + sparse(Iv,Iu,val_JMP,dof_1D,dof_1D);
    
end


%% Transform coeff_mat to wavelet space
vMassV = FMWT_COMP_v * vMassV * FMWT_COMP_v';
GradV  = FMWT_COMP_v * GradV  * FMWT_COMP_v';
FluxV  = FMWT_COMP_v * FluxV  * FMWT_COMP_v';


%% Construct block diagonal for LDG ?
DeltaX = blkdiag(FMWT_COMP_x,FMWT_COMP_x)*...
    DeltaX*...
    blkdiag(FMWT_COMP_x',FMWT_COMP_x');


end