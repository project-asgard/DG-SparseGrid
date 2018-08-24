function [f] = forwardMWT(Lev,Deg,Lmin,Lmax,foo,params)

%% Decompose a 1D function into the multiwavelet basis 

% Get the Forward Multi-Wavelet Transform matrix

FMWT = OperatorTwoScale(Deg,2^Lev);

% Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
% [-1,+1] for performing quadrature.

quad_num = 10;
[quad_x,quad_w] = lgwt(quad_num,-1,1);

% Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% to order k.

p_val = transpose(legendre(quad_x,Deg));

% Get grid spacing.

n=2^(Lev);
h=(Lmax-Lmin)/n;
dof_1D=Deg*n;
f=zeros(dof_1D,1);

for i=0:n-1
    
    % Map quad_x from [-1,+1] to [LMin,LMax] physical domain.
    x = h*(quad_x/2+1/2+i) + Lmin;
    
    % Get the function foo(x) at the quadrature points.
    fxHere = foo(x, params);
    
    this_quad = (quad_w .* fxHere);
       
    % Generate the coefficients for DG basis
    
    % Method 1: Unroll the rows in the GEMM 
    %for thisk=1:Deg
    %    
    %    this_k_legendre = p_val(thisk,:);
    %    f(Deg*i+thisk) = mtimes(this_k_legendre,this_quad);
    %    
    %end
    
    % Method 2: Do the GEMM directly 
    f(i*Deg+1:i*Deg+Deg) = mtimes(p_val,this_quad); 
    
        
end

f = f * h * sqrt(1/h)/2;



% Transfer to multi-DG bases

f = mtimes( FMWT, f );

end