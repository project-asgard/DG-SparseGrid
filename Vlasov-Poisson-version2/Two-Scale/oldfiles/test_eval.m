% Test
clear
close all
% clc
format long e


Deg = 2;
Lstart = 0;
Lend = 1;
Lev = 4;
Lmax = 1;



func = @(x,params)(x);%(sin(pi*x));
% 
% for j_x = 1:Deg
%     for j_y = 1:Deg
%         H1(j_x,j_y) = ((-1)^(j_x+j_y-2)    )*H0(j_x,j_y);
%         G1(j_x,j_y) = ((-1)^(Deg+j_x+j_y-2))*G0(j_x,j_y);
%     end
% end
% 
% quad_num = 20;
% [quad_x,quad_w] = lgwt(quad_num,-1,1);
% 
% % Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% % to order k.
% 
% p_val = transpose(legendre(quad_x,Deg));
% 
% f_coef_MWDG = zeros(Deg*2^Lev,1);
% f_coef_DG = zeros(Deg*2^Lev,1);
% 
% 
% % convert DG coef to MWDG coef
% nx = 2^Lev;
% %--------------------------
% % derive the DG coefficient
% %--------------------------
% for i = 0:nx-1
%     % Map quad_x from [-1,+1] to [Lstart,Lend] physical domain.
%     hx = (Lend-Lstart)/nx;
%     xi_x = hx*(quad_x/2+1/2+i)+Lstart;
%     
%     
%     coef_DG = p_val*(quad_w.*func(xi_x))*hx*sqrt(1/hx)/2;
%     f_coef_DG(Deg*i+1:Deg*(i+1)) = coef_DG;
%     
%     % if plotting, use following
%     % %             Meval_x(2*i+1:2*i+2,Deg*i+1:Deg*(i+1))=sqrt(1/hx)*legendre([-1,1],Deg);
%     % %             x_node(2*i+1:2*i+2)=[i*hx,(i+1)*hx];
%     
% end
% % figure;
% % plot(x_node,Meval_x*f_coef_DG,'b--'); hold on
% % return
% % convert DG to MWDG
% f_tmp_DG=f_coef_DG;
% for loc_lev = Lev:-1:1
%     tmp = [G0 G1]*reshape(f_tmp_DG,2*Deg,2^max(loc_lev-1,0));
%     f_coef_MWDG(Deg*2^max(loc_lev-1,0)+1:Deg*2^loc_lev,1) = tmp(:);
%     
%     tmp = [H0 H1]*reshape(f_tmp_DG,2*Deg,2^max(loc_lev-1,0));
%     f_tmp_DG = tmp(:);
%     
% end
% f_coef_MWDG(1:Deg)=f_tmp_DG;
% 


params.A = 0;
[f_coef_MWDG] = forwardMWT(Lev,Deg,Lstart,Lend,func,params);
% 
% 
% quad_num = 20;
% [quad_x,quad_w] = lgwt(quad_num,-1,1);
% 
% % Get the Legendre basis function evaluated at the Legendre-Gauss nodes up
% % to order k.
% 
% p_val = transpose(legendre(quad_x,Deg));
% 
% f_coef_MWDG = zeros(Deg*2^Lev,1);
% f_coef_DG = zeros(Deg*2^Lev,1);
% 
% 
% % convert DG coef to MWDG coef
% nx = 2^Lev;
% %--------------------------
% % derive the DG coefficient
% %--------------------------
% for i = 0:nx-1
%     % Map quad_x from [-1,+1] to [Lstart,Lend] physical domain.
%     hx = (Lend-Lstart)/nx;
%     xi_x = hx*(quad_x/2+1/2+i)+Lstart;
%     
%     
%     coef_DG = p_val*(quad_w.*func(xi_x))*hx*sqrt(1/hx)/2;
%     f_coef_DG(Deg*i+1:Deg*(i+1)) = coef_DG;
%     
%     % if plotting, use following
% % %             Meval_x(2*i+1:2*i+2,Deg*i+1:Deg*(i+1))=sqrt(1/hx)*legendre([-1,1],Deg);
% % %             x_node(2*i+1:2*i+2)=[i*hx,(i+1)*hx];
%     
% end
% % figure;
% % plot(x_node,Meval_x*f_coef_DG,'b--'); hold on
% % return
% % convert DG to MWDG
% f_tmp_DG=f_coef_DG;
% for loc_lev = Lev:-1:1
%     tmp = [G0 G1]*reshape(f_tmp_DG,2*Deg,2^max(loc_lev-1,0));
%     f_coef_MWDG(Deg*2^max(loc_lev-1,0)+1:Deg*2^loc_lev,1) = tmp(:);
%     
%     tmp = [H0 H1]*reshape(f_tmp_DG,2*Deg,2^max(loc_lev-1,0));
%     f_tmp_DG = tmp(:);
% 
% end
% f_coef_MWDG(1:Deg)=f_tmp_DG;




xx = [Lstart:0.1:Lend]';
% f = zeros(length(xx),1);
% for i = 1:length(xx)
    fval = EvalWavPoint(Lstart,Lend,Lev,Deg,f_coef_MWDG,xx);
%     fval = EvalWavPoint2(Lstart,Lend,Lev,Deg,f_coef_MWDG,xx(i));
%     f(i) = fval;
% end
plot(xx,fval,'r-o',xx,func(xx),'b-<')

figure;
plot(xx,fval-func(xx),'r-o')
[norm(fval-func(xx)) abs(max(fval-func(xx)))]