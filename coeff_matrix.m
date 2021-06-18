function [mat,mat_not_rotated] = coeff_matrix(deg,t,dim,pterm,params, ...
    FMWT_blocks, coeff_level)

% div
%   \int_T u'v dT = \hat{u}v|_{\partial T} - \int_T uv' dT
% Here we shall choose \hat{u} = AVG(u)+JMP(u)/2*cval (cval is given)
% need some test

%% Inputs
% t : time
% dimension : an entry in the pde.dimensions array.
% term_1D : one of the dimension entries in an entry in the pde.terms
% array.

%% Supported boundary condition types (BCL and BCR)
% These are chosen in the pde.dimensions options.
% 'P' == periodic (was 0)
% 'D' == dirichlet (set value of solution) (was 1)
% 'N' == neumann   (set first derivative of solution) (was 2)

%% Available coefficient types (coeff_type)
% 'div'
% 'mass' 
% 'grad' 

%% Note on global vs local Lax-Friedrichs (LF) flux
% We do not (cannot) use local upwinding or LF because selecting
% either the sign of the flow field or the value of the coefficient C could
% be multivalued within the multi-D solution for a single-D coeff_matrix.
% Examples of value LF:: \hat{f} = \hat{A*u} = {{A*u}}+|A|*(1-C)/2*[[u]]
% Denote LF = (1-C)
%   LF = 0 --> Central Flux
%   LF = 1 --> Upwind Flux
%   LF = Max(df/du) -->Lax-Friedrich Flux

%% TODO ...
% * Choice of flux (may require input C)
% * Other BCs are done, but the RHS (with source) needs more work
% * Picking which term type

%%
% pde shortcuts

assert(coeff_level >= dim.lev);

type    = pterm.type;

%%
% dim shortcuts

lev     = coeff_level;

%%
% term shortcuts

pterm.dat   = pterm.dat;
pterm.g       = pterm.g;

%%
% Setup jacobi of variable x and define coeff_mat
N = 2^(lev);
h = (dim.max-dim.min) / N;
dof_1D = deg * N;

%%
%  Get quadrature points and weights.
%  quad_x(:) is quad_num by 1
%  quad_w(:) is quad_num by 1
[quad_x,quad_w] = lgwt(default_quad_number(deg),-1,1);

%%
%  Compute the trace values (values at the left and right of each element for all k)
%  p_L(:) is 1 by deg
%  p_R(:) is 1 by deg
p_L = lin_legendre(-1,deg) * 1/sqrt(h);
p_R = lin_legendre(+1,deg) * 1/sqrt(h);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
p_val  = lin_legendre(quad_x,deg)  * 1/sqrt(h);
Dp_val = lin_dlegendre(quad_x,deg) * 1/sqrt(h) * 2/h;

Jacobi = h/2;

mass = sparse(dof_1D,dof_1D);
div  = sparse(dof_1D,dof_1D);
grad = sparse(dof_1D,dof_1D);

%%
% Convert input dat from wavelet (_W) space to realspace (_R)

if isempty(pterm.dat)
    pterm.dat = ones(dof_1D,1);
end
%dat_R = FMWT' * pterm.dat;
trans_side = 'LT';
dat_R = apply_FMWT_blocks(coeff_level, FMWT_blocks, pterm.dat, trans_side);


%% Loop over all elements in this D
%  Here we construct the 1D coeff_mat in realspace, then transform to
%  wavelet space afterwards.
for i=0:N-1
    
    xL = dim.min + i*h;
    xR = xL + h;
    
    %%
    % Get index ranges for ...
    
    %%
    %  Current element
    c1 = deg*i+1;
    c2 = deg*(i+1);
    c = c1:c2;
    
    %%
    % First element
    first1 = deg*(1-1)+1;
    first2 = deg*((1-1)+1);
    first = first1:first2;
    
    %%
    % Last element
    last1 = deg*(N-1)+1;
    last2 = deg*((N-1)+1);
    last = last1:last2;
    
    %%
    % Map quadrature points from [-1,1] to physical domain of this i element
    
    quad_xi = (((quad_x+1)/2+i)*h+dim.min);
    
    %%
    % Perform volume integral to give deg x deg matrix block
    
    %%
    % Get dat_R at the quadrature points
    
    dat_R_quad = p_val * dat_R(c1:c2);
    
    %%
    % mass matrix u . v
    G1 = pterm.g(quad_xi,params,t,dat_R_quad) .* dim.volume_element(quad_xi,params,t);
    val_mass = p_val' * (G1 .* p_val .* quad_w) * Jacobi;
    
    %%
    % div matrix u . v'
    G1 = pterm.g(quad_xi,params,t,dat_R_quad) .* dim.volume_element(quad_xi,params,t);
    val_div  = -Dp_val'* (G1 .* p_val .* quad_w) * Jacobi;
    
    assert(~isnan(norm(G1)))
    
    Iu = meshgrid( deg*i+1 : deg*(i+1) );
    
    mass = mass + sparse(Iu',Iu,val_mass,dof_1D,dof_1D);
    div = div + sparse(Iu',Iu,val_div,dof_1D,dof_1D);
    
    % this is just here demonstrate a more intuituve way to do this
    % += operation across an aribtray set on indices
    mass_not_rotated = full(mass);
    linear_idx = sub2ind(size(mass_not_rotated),Iu',Iu);
    mass_not_rotated(linear_idx) = val_mass;
    
    div_not_rotated = full(mass);
    linear_idx = sub2ind(size(div_not_rotated),Iu',Iu);
    div_not_rotated(linear_idx) = val_div;
    
    assert(~isnan(sum(mass,'all')))
    assert(~isnan(sum(div,'all')))
    
    if (strcmp(type,'div') || strcmp(type,'grad'))
        
        BCL = pterm.IBCL;
        BCR = pterm.IBCR;
        
        FluxVal = pterm.LF;
        
        %%
        % - <funcCoef*{q},p>
        %----------------------------------------------
        % Numerical Flux is defined as
        % Flux = {{f}} + C/2*[[u]]
        %      = ( f_L + f_R )/2 + FunCoef*( u_R - u_L )/2
        % [[v]] = v_R - v_L
        % F(u) = {{f(u)}} + C/2 * [[u]]
        % {{f(u)}} = (f(u_L) + f(u_R))/2
        % [[u]] = u_R - u_L
        
        FCL = pterm.g(xL,params,t,dat_R_quad) .* dim.volume_element(xL,params,t);
        FCR = pterm.g(xR,params,t,dat_R_quad) .* dim.volume_element(xR,params,t);
        TraVal = [...
            (-p_L)' * FCL * p_R/2 + FluxVal * abs(FCL)/2 * (-p_L)' *   p_R,...
            (-p_L)' * FCL * p_L/2 + FluxVal * abs(FCL)/2 * (-p_L)' * (-p_L),...% xL
            ( p_R)' * FCR * p_R/2 + FluxVal * abs(FCR)/2 * ( p_R)' *   p_R,...
            ( p_R)' * FCR * p_L/2 + FluxVal * abs(FCR)/2 * ( p_R)' * (-p_L),...% xR
            ];
        FCL(isnan(FCL))=0;
        
        
        %%
        % Setup boundary conditions
        
        
        %%
        % If dirichelt
        % u^-_LEFT = g(LEFT)
        % u^+_RIGHT = g(RIGHT)
        
        if strcmp(BCL,'D') %% left dirichlet
            if i==0
                TraVal = [...
                    (-p_L)' * (p_L-p_L),...
                    (-p_L)' * (p_L-p_L),...% xL
                    ( p_R)' * FCR * p_R/2 + FluxVal * abs(FCR)/2 * ( p_R)' *   p_R,...
                    ( p_R)' * FCR * p_L/2 + FluxVal * abs(FCR)/2 * ( p_R)' * (-p_L),...% xR
                    ];
            end
        end
        
        if strcmp(BCR,'D') %% right dirichlet
            if i==N-1
                TraVal = [...
                    (-p_L)' * FCL * p_R/2 + FluxVal * abs(FCL)/2 * (-p_L)' *   p_R,...
                    (-p_L)' * FCL * p_L/2 + FluxVal * abs(FCL)/2 * (-p_L)' * (-p_L),...% xL
                    ( p_R)' * (p_R-p_R),...
                    ( p_R)' * (p_R-p_R),...% xR
                    ];
            end
        end
        
        %%
        % If neumann
        % (divient u)*n = g
        % by splitting div u = q by LDG methods, the B.C is changed to
        % q*n = g (=> q = g for 1D variable)
        % only work for derivatives greater than 1
        
        if strcmp(BCL,'N') %% left neumann
            if i==0
                TraVal = [...
                    (-p_L)' * (p_L-p_L),...,...
                    (-p_L)' * FCL * p_L,...% xL
                    ( p_R)' * FCR * p_R/2 + FluxVal * abs(FCR)/2 * ( p_R)' *   p_R,...
                    ( p_R)' * FCR * p_L/2 + FluxVal * abs(FCR)/2 * ( p_R)' * (-p_L),...% xR
                    ];
                
            end
        end
        
        if strcmp(BCR,'N') %% right neumann
            if i==N-1
                TraVal = [...
                    (-p_L)' * FCL * p_R/2 + FluxVal * abs(FCL)/2 * (-p_L)' *   p_R,...
                    (-p_L)' * FCL * p_L/2 + FluxVal * abs(FCL)/2 * (-p_L)' * (-p_L),...% xL
                    ( p_R)' * FCR * p_R,...
                    ( p_R)' * (p_R-p_R),...% xR
                    ];
                
            end
        end
        
        
        %%
        % Adding trace value to matrix
        
        RowInd = [c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg)];
        ColInd = [ones(deg,1)*(c-deg),ones(deg,1)*c,ones(deg,1)*c,ones(deg,1)*(c+deg)];
        
        if i == 0
            Iu = RowInd(:,deg+1:end);
            Iv = ColInd(:,deg+1:end);
            Val = TraVal(:,deg+1:end);
        elseif i == N - 1
            Iu = RowInd(:,1:3*deg);
            Iv = ColInd(:,1:3*deg);
            Val = TraVal(:,1:3*deg);
        else
            Iu = RowInd;
            Iv = ColInd;
            Val = TraVal;
        end
        
        %%
        % If periodic (Note: the order of this block relative to above matters)
        
        if strcmp(BCL,'P') || strcmp(BCR,'P') %% periodic'
            
            if lev>=1
                if i==0
                    RowInd = [c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg)];
                    ColInd = [ones(deg,1)*last,ones(deg,1)*c,ones(deg,1)*c,ones(deg,1)*(c+deg)];
                end
                if i==N-1
                    RowInd = [c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg)];
                    ColInd = [ones(deg,1)*(c-deg),ones(deg,1)*c,ones(deg,1)*c,ones(deg,1)*first];
                end
            else
                RowInd = [c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg) c'*ones(1,deg)];
                ColInd = [ones(deg,1)*c,ones(deg,1)*c,ones(deg,1)*c,ones(deg,1)*c];
            end
            
            Iu = RowInd;
            Iv = ColInd;
            Val = TraVal;
            
        end
        
        % this is just here demonstrate a more intuituve way to do this
        % += operation across an aribtray set on indices
        linear_idx = sub2ind(size(div_not_rotated),Iu,Iv);
        div_not_rotated(linear_idx) = div_not_rotated(linear_idx) + Val;
        
        div = div + sparse(Iu,Iv,Val,dof_1D,dof_1D);
        assert(~isnan(sum(div,'all')))
        
    end
    
end

grad = -div';

%%
% Store non-transformed matrices for convenince

mass_not_rotated = mass;
div_not_rotated  = div;
grad_not_rotated = grad;

%% Transform coeff_mat to wavelet space
% mat = FMWT * mat * FMWT';

left_notrans = 'LN';
right_trans  = 'RT';

mass = apply_FMWT_blocks(coeff_level, FMWT_blocks, mass, left_notrans);
mass = apply_FMWT_blocks(coeff_level, FMWT_blocks, mass,  right_trans);

div  = apply_FMWT_blocks(coeff_level, FMWT_blocks,  div, left_notrans);
div  = apply_FMWT_blocks(coeff_level, FMWT_blocks,  div,  right_trans);

grad = apply_FMWT_blocks(coeff_level, FMWT_blocks, grad, left_notrans);
grad = apply_FMWT_blocks(coeff_level, FMWT_blocks, grad,  right_trans);

assert(~isnan(sum(mass,'all')))
assert(~isnan(sum( div,'all')))
assert(~isnan(sum(grad,'all')))

if strcmp(type, 'div')
    mat  = div;
    mat_not_rotated = div_not_rotated;
end
if strcmp(type,'mass')
    mat  = mass;
    mat_not_rotated = mass_not_rotated;
end
if strcmp(type,'grad')
    mat  = grad;
    mat_not_rotated = grad_not_rotated;
end

end
