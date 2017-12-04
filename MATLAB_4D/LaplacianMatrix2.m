% function [M_stiff,b,xnode,Meval,uu]=LaplacianMatrix2(n,k)
function [L_MWDG,RHS_MWDG,Meval]=LaplacianMatrix2(n,k)
% PDE dependent Stiffness Matrix
% working version with two-scale relation

% clc
% clear all
% close all
% 
% n=5; k=2;
%--DG parameters
h=1/2^n;
theta=1;
sigma=10*k;
quad_num=50;
%---------------

dof = 2^n*k;

% compute the trace values
Dp_1 = dlegendre(-1,k);
Dp_2 = dlegendre( 1,k);
p_1 = legendre(-1,k);
p_2 = legendre(1,k);

[quad_x,quad_w]=lgwt(quad_num,-1,1);

p_val = legendre(quad_x,k);
Dp_val = dlegendre(quad_x,k);

 
% % <{grad u},[v]>
% [-p_1'*Dp_2/2,-p_1'*Dp_1/2+p_2'*Dp_2,p_2'*Dp_1/2] 
% % <{grad v},[u]>
% [Dp_1'*p_2/2,-Dp_1'*p_1/2+Dp_2'*p_2/2,-Dp_2'*p_1/2]
% % <[u],[v]>
% [-p_1'*p_2,p_1'*p_1+p_2'*p_2,-p_2'*p_1]

dof_1D=k*2^n;
M_stiff=sparse(dof_1D,dof_1D);
b=sparse(dof_1D,1);

% generate 1D stiffness matrix for DG 
for LL=0:2^n-1
    
    if LL==0
        Iu=[meshgrid(k*LL+1:k*(LL+1)),meshgrid(k*(LL+1)+1:k*(LL+2))];
        Iv=[meshgrid(k*LL+1:k*(LL+1))',meshgrid(k*LL+1:k*(LL+1))'];
        val=2^(2*n+1)*[Dp_val'*(quad_w.*Dp_val)+...
            -(-p_1'*Dp_1+p_2'*Dp_2/2)-(-Dp_1'*p_1+Dp_2'*p_2/2)+sigma/2*(p_1'*p_1+p_2'*p_2),...
            -(p_2'*Dp_1/2)-(-Dp_2'*p_1/2)+sigma/2*(-p_2'*p_1)...
            ];
        
    elseif LL==2^n-1
        Iu=[meshgrid(k*(LL-1)+1:k*(LL)),meshgrid(k*LL+1:k*(LL+1))];
        Iv=[meshgrid(k*LL+1:k*(LL+1))',meshgrid(k*LL+1:k*(LL+1))'];
        val=2^(2*n+1)*[-(-p_1'*Dp_2/2)-(Dp_1'*p_2/2)+sigma/2*(-p_1'*p_2),...
            Dp_val'*(quad_w.*Dp_val)+...
            -(-p_1'*Dp_1/2+p_2'*Dp_2)-(-Dp_1'*p_1/2+Dp_2'*p_2)+sigma/2*(p_1'*p_1+p_2'*p_2)...
            ];
        
    else
        Iu=[meshgrid(k*(LL-1)+1:k*(LL)),meshgrid(k*LL+1:k*(LL+1)),meshgrid(k*(LL+1)+1:k*(LL+2))];
        Iv=[meshgrid(k*LL+1:k*(LL+1))', meshgrid(k*LL+1:k*(LL+1))',meshgrid(k*LL+1:k*(LL+1))'];
        val=2^(2*n+1)*[-(-p_1'*Dp_2/2)-(Dp_1'*p_2/2)+sigma/2*(-p_1'*p_2),...
            Dp_val'*(quad_w.*Dp_val)+...
            -(-p_1'*Dp_1/2+p_2'*Dp_2/2)-(-Dp_1'*p_1/2+Dp_2'*p_2/2)+sigma/2*(p_1'*p_1+p_2'*p_2),...
             -(p_2'*Dp_1/2)-(-Dp_2'*p_1/2)+sigma/2*(-p_2'*p_1);
            ];        
    end
    
    M_stiff=M_stiff+sparse(Iu,Iv,val,dof_1D,dof_1D);
    
    % RHS term
     ff=rhs( quad_x*2^(-n-1)+2^(-n)*(LL+0.5) );
     Iu=[k*LL+1:k*(LL+1)];
     b(Iu)=p_val'*(quad_w.*ff)*2^(-n)/2*2^(n/2);
     
     
     Meval(2*LL+1:2*LL+2,Iu)=2^(n/2)*legendre([-1,1],k);

end

% full(M_stiff)
% 
% return

sol_1D=M_stiff\b*pi^2;

% % max(sol_1D-b)
% % subplot(1,2,1)
% % plot(Meval*sol_1D)
% % hold on
% % plot(Meval*b,'ro')

% Transfer by Two-Scale Relationship
load(['two_scale_rel_',num2str(k),'.mat'])
Np=2^n;
H0(find(abs(H0)<1e-10))=0;
G0(find(abs(G0)<1e-10))=0;
%------------------------------
% Set-up Diffusion operator %%
%------------------------------

H1 = zeros(k);
G1 = zeros(k);

for j_x = 1:k
    for j_y = 1:k
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2))*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(k+j_x+j_y-2))*G0(j_x,j_y);
    end
end

FMWT = zeros(k*Np);
iFMWT = zeros(k*Np);

for j=1:Np/2
    % The reverse order from Lin
    FMWT(k*(j-1)+1:k*j,2*k*(j-1)+1:2*k*j)=[H0 H1];
    FMWT(k*(j+Np/2-1)+1:k*(j+Np/2),2*k*(j-1)+1:2*k*j) = [G0 G1];
end
iFMWT=FMWT';

sp = [];
FMWT_COMP = eye(k*Np);
for j=1:n
    cFMWT = FMWT;
    % Reverse the index in matrix from Lin
    if j>1
        cFMWT = zeros(k*Np);
        cn = 2^(n-j+1)*k;
        cnr=Np*k-cn;
        cFMWT(cn+1:k*Np,cn+1:k*Np)=eye(Np*k-cn);
        cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
        cFMWT(cn/2+1:cn,1:cn)=FMWT(k*Np/2+1:k*Np/2+cn/2,1:cn);
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


