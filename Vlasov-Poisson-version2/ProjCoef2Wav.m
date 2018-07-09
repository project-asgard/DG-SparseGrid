function FF = ProjCoef2Wav(LevX,Deg,Lmax,FMWT_COMP_x,func)
% projection from a given function to the Multiwavelet functions
% Pro(func) = sum F_i*phi_i
% FF = [F_i]

% Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
% [-1,+1] for performing quadrature.

quad_num = 10;
[quad_x,quad_w] = lgwt(quad_num,-1,1);

% Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% to order k.

p_val = transpose(legendre(quad_x,Deg));


% Get grid spacing for both  x and v.

nx=2^(LevX);
hx=Lmax/nx;
dof_1D=Deg*nx;
f_x=zeros(dof_1D,1);



for i=0:nx-1
    
    % Map quad_x from [-1,+1] to [0,LMax] physical domain.
    xi_x = hx*(quad_x/2+1/2+i);
    
    Meval_x(2*i+1:2*i+2,Deg*i+1:Deg*(i+1))=sqrt(1/hx)*legendre([-1,1],Deg);
    x_node(2*i+1:2*i+2)=[i*hx,(i+1)*hx];
    
%      f_x(Deg*i+1:Deg*(i+1))=p_val*(quad_w.*func(xi_x))*hx*sqrt(1/hx)/2;
     
    % Get the f(x) at the quadrature points.
    fxHere = func(xi_x);
       
    % Generate the coefficients for DG basis    
    for thisk=1:Deg
        
        this_k_legendre = p_val(thisk,:);
        this_quad = (quad_w .* fxHere);
        f_x(Deg*i+thisk) = mtimes(this_k_legendre,this_quad) * hx * sqrt(1/hx)/2;
        
    end
    
end

plot(x_node,Meval_x*f_x,'r-o',x_node,func(x_node),'b--')
% Transfer to multi-DG bases

FF = mtimes( FMWT_COMP_x, f_x );
figure
plot(x_node,Meval_x*FMWT_COMP_x'*FF,'r-o',x_node,func(x_node),'b--')

end