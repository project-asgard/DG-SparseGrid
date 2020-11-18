function f_coef_MWDG = proj_coef_to_wav2(Lev,Deg,Lstart,Lend,func)
%f_coef_MWDG = ProjCoef2Wav_v2(Lev,Deg,Lmax,func)
% projection from a given function to the Multiwavelet functions
% Pro(func) = sum F_i*phi_i
% FF = [F_i]

% Get the Legendre-Gauss nodes (quad_x) and weights (quad_w) on the domain
% [-1,+1] for performing quadrature.
load(['two_scale_rel_',num2str(Deg),'.mat'])

for j_x = 1:Deg
    for j_y = 1:Deg
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)    )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(Deg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

quad_num = 10;
[quad_x,quad_w] = lgwt(quad_num,-1,1);

% Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% to order k.

p_val = transpose(legendre(quad_x,Deg));

f_coef_MWDG = zeros(Deg*2^Lev,1);
f_coef_DG = zeros(Deg*2^Lev,1);
f = zeros(Deg*2^Lev,1);

% convert DG coef to MWDG coef
nx = 2^Lev;
%--------------------------
% derive the DG coefficient
%--------------------------
for i = 0:nx-1
    % Map quad_x from [-1,+1] to [Lstart,Lend] physical domain.
    hx = (Lend-Lstart)/nx;
    xi_x = hx*(quad_x/2+1/2+i)+Lstart;
    
    this_quad = quad_w.*func(xi_x);
    coef_DG = (p_val * this_quad) * hx*sqrt(1/hx)/2;
    f_coef_DG(Deg*i+1:Deg*(i+1)) = coef_DG;
    
    fxHere = func(xi_x);
    for thisk=1:Deg
        
        this_k_legendre = p_val(thisk,:);
        %this_quad = (quad_w .* fxHere);
        f(Deg*i+thisk) = mtimes(this_k_legendre,this_quad) * hx * sqrt(1/hx)/2;
        %f(Deg*i+thisk) = this_k_legendre * this_quad;%* hx * sqrt(1/hx)/2;
        %f(Deg*i+thisk) = sum(this_k_legendre .* this_quad');

        
    end
    
    % if plotting, use following
    % %         Meval_x(2*i+1:2*i+2,Deg*i+1:Deg*(i+1))=sqrt(1/hx)*legendre([-1,1],Deg);
    % %         x_node(2*i+1:2*i+2)=[i*hx,(i+1)*hx];
    
end

% convert DG to MWDG
f_tmp_DG=f_coef_DG;
for loc_lev = Lev:-1:1
    tmp = [G0 G1]*reshape(f_tmp_DG,2*Deg,2^max(loc_lev-1,0));
    f_coef_MWDG(Deg*2^max(loc_lev-1,0)+1:Deg*2^loc_lev,1) = tmp(:);
    
    tmp = [H0 H1]*reshape(f_tmp_DG,2*Deg,2^max(loc_lev-1,0));
    f_tmp_DG = tmp(:);

end
f_coef_MWDG(1:Deg)=f_tmp_DG;


% convert MWDG coef to DG coef
% from MWDG to DG
tmp_DG = f_coef_MWDG(1:Deg);
for loc_lev=1:Lev
    
    tmp = [H0' G0';H1' G1']*[reshape(tmp_DG,Deg,2^(loc_lev-1));...
        reshape(f_coef_MWDG(Deg*2^(loc_lev-1)+1:Deg*2^loc_lev),Deg,2^(loc_lev-1))];
    tmp_DG = tmp(:);
end

% % % by checking
% % FMWT_COMP = OperatorTwoScale(Deg,2^Lev);
% % subplot(1,2,1)
% % plot(tmp_DG-f_coef_DG)
% % subplot(1,2,2)
% % plot(f_coef_MWDG-FMWT_COMP*f_coef_DG)

end
