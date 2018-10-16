function [L_MWDG,RHS_MWDG,Meval,coef_MW]=LaplacianMatrix2(n,Deg,k)
% PDE dependent Stiffness Matrix
% worDeging version with two-scale relation
% k - Wave Number
%--DG parameters
h=1/2^n;
theta=1;
sigma=10*Deg^2*k^2;
quad_num=50;
%---------------

dof = 2^n*Deg;

% compute the trace values
Dp_1 = dlegendre(-1,Deg);
Dp_2 = dlegendre( 1,Deg);
p_1 = legendre(-1,Deg);
p_2 = legendre(1,Deg);

[quad_x,quad_w]=lgwt(quad_num,-1,1);

p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);


rhs = @(x,k)((pi^2-k^2)*sin(pi*x));
exactu = @(x,k)(sin(pi*x));
bc_g = @(x,k,n)(pi*cos(pi*x)*n+1i*k*sin(pi*x));

% rhs = @(x,k)((k^2*pi^2-k^2)*sin(pi*k*x));
% exactu = @(x,k)(sin(pi*k*x));
% bc_g = @(x,k,n)(k*pi*cos(k*pi*x)*n+1i*k*sin(pi*k*x));

eps = k^2;
eta = sqrt(k^2+1i*eps);
eta2 = sqrt(k^2);%eta;%

% % <{grad u},[v]>
% [-p_1'*Dp_2/2,-p_1'*Dp_1/2+p_2'*Dp_2,p_2'*Dp_1/2] 
% % <{grad v},[u]>
% [Dp_1'*p_2/2,-Dp_1'*p_1/2+Dp_2'*p_2/2,-Dp_2'*p_1/2]
% % <[u],[v]>
% [-p_1'*p_2,p_1'*p_1+p_2'*p_2,-p_2'*p_1]

dof_1D=Deg*2^n;
M_stiff=sparse(dof_1D,dof_1D);
Aeps = sparse(dof_1D,dof_1D);

b=sparse(dof_1D,1);
coef_DG=sparse(dof_1D,1);
% generate 1D stiffness matrix for DG 
for LL=0:2^n-1
    
        
%        val=2^(2*n+1)*[-(-p_1'*Dp_2/2)-(Dp_1'*p_2/2)+sigma/2*(-p_1'*p_2),... % previous
%             Dp_val'*(quad_w.*Dp_val)+...
%             -(-p_1'*Dp_1/2+p_2'*Dp_2/2)-(-Dp_1'*p_1/2+Dp_2'*p_2/2)+sigma/2*(p_1'*p_1+p_2'*p_2),...
%              -(p_2'*Dp_1/2)-(-Dp_2'*p_1/2)+sigma/2*(-p_2'*p_1); % later cell
%             ]; 
       xi = quad_x*2^(-n-1)+2^(-n)*(LL+0.5);
  
       tmp_val=2^(2*n+1)*[...
           -(-p_1'*Dp_2/2)-(Dp_1'*p_2/2)+sigma/2*(-p_1'*p_2),... % previous
           ...
            Dp_val'*(quad_w.*Dp_val)+...
            -k^2*eye(Deg)/(2^(2*n+1))+...
            -(-p_1'*Dp_1/2+p_2'*Dp_2/2)-(-Dp_1'*p_1/2+Dp_2'*p_2/2)+sigma/2*(p_1'*p_1+p_2'*p_2),...
           ... 
            -(p_2'*Dp_1/2)-(-Dp_2'*p_1/2)+sigma/2*(-p_2'*p_1); % later cell
            ]; 
        
    if LL==0 % Initial Cell
        Iu=[meshgrid(Deg*LL+1:Deg*(LL+1)),meshgrid(Deg*(LL+1)+1:Deg*(LL+2))];
        Iv=[meshgrid(Deg*LL+1:Deg*(LL+1))',meshgrid(Deg*LL+1:Deg*(LL+1))'];
        val = tmp_val (:,Deg+1:end);
        
        val_bc = 1i*k*(-p_1)'*(p_1)/2*2^(2*n+1);
        
        val(:,1:Deg) = val(:,1:Deg) + val_bc;
        
        val_bc = -p_1' * bc_g(0,k,-1)*2^(-n)/2*2^(n/2);
    elseif LL==2^n-1 % End Cell
        Iu=[meshgrid(Deg*(LL-1)+1:Deg*(LL)),meshgrid(Deg*LL+1:Deg*(LL+1))];
        Iv=[meshgrid(Deg*LL+1:Deg*(LL+1))',meshgrid(Deg*LL+1:Deg*(LL+1))'];
        val = tmp_val(:,1:2*Deg);
        val_bc = 1i*k*(-p_1)'*(p_1)/2*2^(2*n+1);
        
        val(:,Deg+1:2*Deg) = val(:,Deg+1:2*Deg) + val_bc;
        val_bc = p_2' * bc_g(1,k,1)*2^(-n)/2*2^(n/2);
    else
        Iu=[meshgrid(Deg*(LL-1)+1:Deg*(LL)),meshgrid(Deg*LL+1:Deg*(LL+1)),meshgrid(Deg*(LL+1)+1:Deg*(LL+2))];
        Iv=[meshgrid(Deg*LL+1:Deg*(LL+1))', meshgrid(Deg*LL+1:Deg*(LL+1))',meshgrid(Deg*LL+1:Deg*(LL+1))']; 
        val = tmp_val;
        val_bc = zeros(Deg,1);
    end
    
    M_stiff=M_stiff+sparse(Iu,Iv,val,dof_1D,dof_1D);
    
    % RHS term
    
    ff=rhs( xi,k );
    Iu=[Deg*LL+1:Deg*(LL+1)];
    b(Iu)=p_val'*(quad_w.*ff)*2^(-n)/2*2^(n/2)+val_bc;
     
     
     
     Meval(2*LL+1:2*LL+2,Iu)=2^(n/2)*legendre([-1,1],Deg);
     uu=exactu( xi,k );
     coef_DG(Iu)=p_val'*(quad_w.*uu)*2^(-n)/2*2^(n/2);
     
   % for Aeps  
     tmp_val=2^(2*n+1)*[...
           -(-p_1'*Dp_2/2)-(Dp_1'*p_2/2)+sigma/2*(-p_1'*p_2),... % previous
           ...
            Dp_val'*(quad_w.*Dp_val)+...
            -eta^2*eye(Deg)/(2^(2*n+1))+...
            -(-p_1'*Dp_1/2+p_2'*Dp_2/2)-(-Dp_1'*p_1/2+Dp_2'*p_2/2)+sigma/2*(p_1'*p_1+p_2'*p_2),...
           ... 
            -(p_2'*Dp_1/2)-(-Dp_2'*p_1/2)+sigma/2*(-p_2'*p_1); % later cell
            ]; 
        
    if LL==0 % Initial Cell
        Iu=[meshgrid(Deg*LL+1:Deg*(LL+1)),meshgrid(Deg*(LL+1)+1:Deg*(LL+2))];
        Iv=[meshgrid(Deg*LL+1:Deg*(LL+1))',meshgrid(Deg*LL+1:Deg*(LL+1))'];
        val = tmp_val (:,Deg+1:end);
        
        val_bc = 1i*eta2*(-p_1)'*(p_1)/2*2^(2*n+1);
        
        val(:,1:Deg) = val(:,1:Deg) + val_bc;
        
%         val_bc = -p_1' * bc_g2(0,k,-1)*2^(-n)/2*2^(n/2);
    elseif LL==2^n-1 % End Cell
        Iu=[meshgrid(Deg*(LL-1)+1:Deg*(LL)),meshgrid(Deg*LL+1:Deg*(LL+1))];
        Iv=[meshgrid(Deg*LL+1:Deg*(LL+1))',meshgrid(Deg*LL+1:Deg*(LL+1))'];
        val = tmp_val(:,1:2*Deg);
        val_bc = 1i*eta2*(-p_1)'*(p_1)/2*2^(2*n+1);
        
        val(:,Deg+1:2*Deg) = val(:,Deg+1:2*Deg) + val_bc;
%         val_bc = p_2' * bc_g(1,k,1)*2^(-n)/2*2^(n/2);
    else
        Iu=[meshgrid(Deg*(LL-1)+1:Deg*(LL)),meshgrid(Deg*LL+1:Deg*(LL+1)),meshgrid(Deg*(LL+1)+1:Deg*(LL+2))];
        Iv=[meshgrid(Deg*LL+1:Deg*(LL+1))', meshgrid(Deg*LL+1:Deg*(LL+1))',meshgrid(Deg*LL+1:Deg*(LL+1))']; 
        val = tmp_val;
%         val_bc = zeros(Deg,1);
    end
    
    Aeps = Aeps+sparse(Iu,Iv,val,dof_1D,dof_1D);
    

end

% full(M_stiff)
% 
% return



MaxIter = min(1000,size(M_stiff,1));
Tol = 1e-10;
Max_NUM = min(100,size(M_stiff,1));

[sol_GMRES,fl0,rr0,it0,rv_RMRES]= gmres(M_stiff,b,Max_NUM,Tol,MaxIter,Aeps);
subplot(1,2,2)
semilogy(0:length(rv_RMRES)-1,rv_RMRES/norm(ff),'-o','linewidth',2,'markersize',8);
legend({'GMRES'},'fontsize',20)
xlabel('Iteration number');
ylabel('Relative residual');
[max(abs(sol_GMRES-coef_DG)) norm(sol_GMRES-coef_DG)]
it0(2)

sol_1D=M_stiff\b;%*pi^2;

max(sol_1D-coef_DG)
subplot(1,2,1)
plot(Meval*real(sol_1D),'b+')
hold on
% plot(Meval*imag(sol_1D),'r+')
plot(Meval*coef_DG,'ro')

[condest(M_stiff),condest(Aeps\M_stiff)]

% Transfer by Two-Scale Relationship
load(['two_scale_rel_',num2str(Deg),'.mat'])
Np=2^n;
H0(find(abs(H0)<1e-10))=0;
G0(find(abs(G0)<1e-10))=0;
%------------------------------
% Set-up Diffusion operator %%
%------------------------------

H1 = zeros(Deg);
G1 = zeros(Deg);

for j_x = 1:Deg
    for j_y = 1:Deg
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2))*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(Deg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

FMWT = zeros(Deg*Np);
iFMWT = zeros(Deg*Np);

for j=1:Np/2
    % The reverse order from Lin
    FMWT(Deg*(j-1)+1:Deg*j,2*Deg*(j-1)+1:2*Deg*j)=[H0 H1];
    FMWT(Deg*(j+Np/2-1)+1:Deg*(j+Np/2),2*Deg*(j-1)+1:2*Deg*j) = [G0 G1];
end
iFMWT=FMWT';

sp = [];
FMWT_COMP = eye(Deg*Np);
for j=1:n
    cFMWT = FMWT;
    % Reverse the index in matrix from Lin
    if j>1
        cFMWT = zeros(Deg*Np);
        cn = 2^(n-j+1)*Deg;
        cnr=Np*Deg-cn;
        cFMWT(cn+1:Deg*Np,cn+1:Deg*Np)=eye(Np*Deg-cn);
        cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
        cFMWT(cn/2+1:cn,1:cn)=FMWT(Deg*Np/2+1:Deg*Np/2+cn/2,1:cn);
    end
%     j
%     cFMWT

    FMWT_COMP = cFMWT*FMWT_COMP;
%     FMWT_COMP
end

iFMWT_COMP=FMWT_COMP';

L_MWDG = FMWT_COMP*M_stiff*iFMWT_COMP;
RHS_MWDG=FMWT_COMP*b;%*pi^2;

% % Truncate the entry in the matrix of MWDG
epsilon_MWDG=0;
L_MWDG(find(abs(L_MWDG)<epsilon_MWDG))=0;

% Meval=FMWT_COMP*Meval;

Meval=Meval*iFMWT_COMP;
% 
coef_MW=FMWT_COMP*coef_DG;



% figure
% plot(Meval*(coef_DG))
% hold on
% plot(Meval*b,'ro')
return
% Compute the solution from MWDG
sol_MW=(L_MWDG\RHS_MWDG);

% full([sol_1D iFMWT_COMP*sol_MW])

subplot(1,2,2)
plot(Meval*(iFMWT_COMP*sol_MW))
hold on
plot(Meval*b,'ro')

figure;
spy(L_MWDG)
