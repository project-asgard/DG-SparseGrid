function FF = ProjCoef2Wav(LevX,Deg,Lmax,func)
% projection from a given function to the Multiwavelet functions
% Pro(func) = sum F_i*phi_i
% FF = [F_i]

% Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
% [-1,+1] for performing quadrature.
load(['two_scale_rel_',num2str(Deg),'.mat'])

for j_x = 1:Deg
    for j_y = 1:Deg
        G1(j_x,j_y) = ((-1)^(Deg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

FMWT_COMP_x = OperatorTwoScale(Deg,2^LevX);

quad_num = 10;
[quad_x,quad_w] = lgwt(quad_num,-1,1);

% Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% to order k.

p_val = transpose(legendre(quad_x,Deg));

for loc_lev = 0:LevX
    nx = 2^loc_lev;
    %--------------------------
    % derive the DG coefficient
    %--------------------------
    for i = 0:nx-1 
        % Map quad_x from [-1,+1] to [0,LMax] physical domain.
        hx = Lmax/nx;
        xi_x = hx*(quad_x/2+1/2+i);
     
        coef_DG = p_val*(quad_w.*func(xi_x))*hx*sqrt(1/hx)/2;
        
        f_coef_DG(Deg*i+1:Deg*(i+1)) = coef_DG;
        
        Meval_x(2*i+1:2*i+2,Deg*i+1:Deg*(i+1))=sqrt(1/hx)*legendre([-1,1],Deg);
        x_node(2*i+1:2*i+2)=[i*hx,(i+1)*hx];
       
    end
    %----------------------------
    % derive the MWDG coefficient
    %----------------------------
        if loc_lev == 0
            f_coef_MWDG(1:Deg) = f_coef_DG;
        else
            
            for i = 0:2^max(0,loc_lev-1)-1
                index_DG1 = Deg*(2*i)+[1:Deg];
                coef1 = f_coef_DG(index_DG1);
                index_DG2 = Deg*(2*i+1)+[1:Deg];
                coef2 = f_coef_DG(index_DG2);

                index_MWDG = Deg*2^(loc_lev-1)+[Deg*i+1:Deg*(i+1)];
                f_coef_MWDG(index_MWDG) = G0*coef1'+G1*coef2';        
            end
        end
    
end


size(Meval_x)
size(f_coef_DG)
plot(x_node,Meval_x*f_coef_DG','r-o',x_node,func(x_node),'b-+')

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

% plot(x_node,Meval_x*f_x,'r-o',x_node,func(x_node),'b--')
% Transfer to multi-DG bases

FF = mtimes( FMWT_COMP_x, f_x );
[FF f_coef_MWDG']
figure
plot(FF-f_coef_MWDG')
% figure
% plot(x_node,Meval_x*FMWT_COMP_x'*FF,'r-o',x_node,func(x_node),'b--')


end