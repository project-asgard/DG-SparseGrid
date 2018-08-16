function [f,nodes,M] = forwardMWT(Lev,Deg,Lmin,Lmax,foo,params,nQuad)

%% Decompose a 1D function into the multiwavelet basis

% Get the Forward Multi-Wavelet Transform matrix

FMWT = OperatorTwoScale(Deg,2^Lev);

% Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
% [-1,+1] for performing quadrature.

if ~exist('nQuad','var') || isempty(nQuad)
    quad_num = 10;
else
    quad_num = nQuad;
end

[quad_x,quad_w] = lgwt(quad_num,-1,1);

% Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% to order k.

p_val = transpose(legendre(quad_x,Deg));

% Get grid spacing.

n=2^(Lev);
h=(Lmax-Lmin)/n;
dof_1D=Deg*n;
nodes=zeros(dof_1D,1);
f=zeros(dof_1D,1);
M=zeros(dof_1D,dof_1D);

for i=0:n-1
    
    % Map quad_x from [-1,+1] to [LMin,LMax] physical domain.
    x = h*(quad_x/2+1/2+i) + Lmin;
    
    % Get the function foo(x) at the quadrature points.
    fxHere = foo(x, params);
    
    this_quad = (quad_w .* fxHere);
    
    % Generate the coefficients for DG basis
    
    % Do the GEMM directly
    f(i*Deg+1:i*Deg+Deg) = mtimes(p_val,this_quad);
    
    
    if(quad_num == Deg) % Why is this the case only for the return transform?
        % Store the locations of the nodes
        nodes(i*Deg+1:i*Deg+Deg) = x;
        
        % Store the Lengende functions for use in transforming back
        M(i*Deg+1:i*Deg+Deg,i*Deg+1:i*Deg+Deg) = p_val;
    end
end

f = f * h * sqrt(1/h)/2;
M = M * sqrt(1/h);

% Transfer to multi-DG bases

f = mtimes( FMWT, f );
M = mtimes( FMWT, M );
M = M';

end