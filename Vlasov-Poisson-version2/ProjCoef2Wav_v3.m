function f_coef_MWDG = ProjCoef2Wav_v3(Lev,Deg,Lstart,Lend,func)
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

quad_num = 20;
[quad_x,quad_w] = lgwt(quad_num,-1,1);

% Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% to order k.

p_val = transpose(legendre(quad_x,Deg));

f_coef_MWDG = zeros(Deg*2^Lev,1);
for loc_lev = 0:Lev
    nx = 2^loc_lev;
    %--------------------------
    % derive the DG coefficient
    %--------------------------
    f_coef_DG = zeros(Deg*nx,1);
    for i = 0:nx-1
        % Map quad_x from [-1,+1] to [Lstart,Lend] physical domain.
        %hx = Lmax/nx;
        %         xi_x = hx*(quad_x/2+1/2+i);
        
        hx = (Lend-Lstart)/nx;
        xi_x = hx*(quad_x/2+1/2+i)+Lstart;
        
        
        coef_DG = p_val*(quad_w.*func(xi_x))*hx*sqrt(1/hx)/2;
        
        f_coef_DG(Deg*i+1:Deg*(i+1)) = coef_DG;
        
        % if plotting, use following
        % %         Meval_x(2*i+1:2*i+2,Deg*i+1:Deg*(i+1))=sqrt(1/hx)*legendre([-1,1],Deg);
        % %         x_node(2*i+1:2*i+2)=[i*hx,(i+1)*hx];
        
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
            f_coef_MWDG(index_MWDG) = G0*coef1+G1*coef2;
        end
    end
    
end

% compute the finest lev for DG coeff
% convert DG coef to MWDG coef
nx = 2^Lev;
%--------------------------
% derive the DG coefficient
%--------------------------
f_coef_DG = zeros(Deg*nx,1);
for i = 0:nx-1
    % Map quad_x from [-1,+1] to [Lstart,Lend] physical domain.
    
    hx = (Lend-Lstart)/nx;
    xi_x = hx*(quad_x/2+1/2+i)+Lstart;
    
    
    coef_DG = p_val*(quad_w.*func(xi_x))*hx*sqrt(1/hx)/2;
    
    f_coef_DG(Deg*i+1:Deg*(i+1)) = coef_DG;
    
    % if plotting, use following
    % %         Meval_x(2*i+1:2*i+2,Deg*i+1:Deg*(i+1))=sqrt(1/hx)*legendre([-1,1],Deg);
    % %         x_node(2*i+1:2*i+2)=[i*hx,(i+1)*hx];
    
end

% % % for loc_lev = Lev:-1:1
% % %     tmp = [G0 G1]*reshape(f_coef_DG,2*Deg,2^max(loc_lev-1,0));
% % %     f_tmp_MWDG(Deg*2^max(loc_lev-1,0)+1:Deg*2^loc_lev,1) = tmp(:);
% % %     
% % %     tmp = [H0 H1]*reshape(f_coef_DG,2*Deg,2^max(loc_lev-1,0));
% % %     f_tmp_DG = tmp(:);
% % %     
% % %     f_coef_DG = f_tmp_DG;
% % % end
% % % f_tmp_MWDG(1:Deg)=f_tmp_DG;
% % % 
% % % % [f_tmp_MWDG,f_coef_MWDG]
% % % 
% % % plot(f_tmp_MWDG-f_coef_MWDG)


% % size(Meval_x)
% % size(f_coef_DG)
% % plot(x_node,Meval_x*f_coef_DG,'r-o',x_node,func(x_node),'b-+')

% convert MWDG coef to DG coef

% from MWDG to DG
tmp_DG = f_coef_MWDG(1:Deg);
for loc_lev=1:Lev
    
    tmp = [H0' G0';H1' G1']*[reshape(tmp_DG,Deg,2^(loc_lev-1));...
        reshape(f_coef_MWDG(Deg*2^(loc_lev-1)+1:Deg*2^loc_lev),Deg,2^(loc_lev-1))];
    tmp_DG = tmp(:);
end


plot(tmp_DG-f_coef_DG)

% % coef_DG_lowlev(1:Deg,1) = f_coef_MWDG(1:Deg);
% % coef_DG(1:Deg,1)=coef_DG_lowlev;
% % for loc_lev = 1:Lev
% %     for i = 0:2^max(0,loc_lev-1)-1
% %
% %         index_DG = Deg*i+[1:Deg];
% %          index_MWDG1 = Deg*2^(loc_lev-1)+i*Deg+[1:Deg];
% %
% %         index_DG1 = Deg*2*i+[1:Deg];
% %         index_DG2 = Deg*2*i+Deg+[1:Deg];
% %
% %         coef1 = coef_DG_lowlev(index_DG);
% %         coef2 = f_coef_MWDG(index_MWDG1);
% %
% %         coef = H0'*coef1+G0'*coef2;
% %         coef_DG(index_DG1) = coef;
% %
% %         coef = H1'*coef1+G1'*coef2;
% %         coef_DG(index_DG2) = coef;
% %     end
% %     coef_DG_lowlev=coef_DG;
% % end

% % [coef_DG f_coef_DG]
% % plot(coef_DG- f_coef_DG)
end