function [M_stiff,b,xnode,Meval,uu]=LaplacianMatrix(n,k)
% PDE dependent Stiffness Matrix


%--DG parameters
h=1/2^n;
theta=1;
sigma=100*h^(-1);
quad_num=10;
%---------------

dof = 2^n*k;

s1=1;
s2=1;
h=1/2^n;

M_stiff=sparse(dof,dof);
M_stiff_diag=sparse(dof,dof);
M_mass=sparse(dof,dof);
M_mass_diag=sparse(dof,dof);
b=zeros(dof,1);

xnode=[0:h:1];
xnode=[xnode(1:end-1);xnode(2:end)];
xnode=xnode(:);
Meval=sparse(size(xnode,1),dof);

for n1 = 0:n
    for l1=0:2^max(0,n1-1)-1
        % consider the other basis phi_2 based on cell I^{n2}_{l2}
        % if the considered basis based on the subset of cell: I^{n2}_{l2}\in
        % I^{n1}_{l1}
        % generate the Gaussian quadrature on I^{n2}_{l2}=2^(-n2)*(l2,l2+1)
        for  n2=n1:n
            
            ratio=2^max(0,n2-n1);
            
            for l2=max(0,ratio*l1-1):min(min(ratio*(l1+1)-1,2^max(0,n2-1)-1)+1,2^max(0,n2-1)-1)
                
                % generate the integration cell
                x1=2^(-n2+1)*l2;
                if n2==0
                    x2=1;
                else
                    x2=2^(-n2+1)*(l2+1);
                end
                
                
                xm=(x1+x2)/2;
                % divide the whole to two different sub-cells
                
                [quad_x,quad_w]=lgwt(quad_num,x1,(x1+x2)/2);
                [quad_x2,quad_w2]=lgwt(quad_num,(x1+x2)/2,x2);
                
                quad_x=[quad_x;quad_x2];
                quad_w=[quad_w;quad_w2];
                
                if n1 ==0
                    % generate the basis phi_1
                    
                    phi_1  = legendre(2*quad_x-1,k);
                    dphi_1 = 2*dlegendre(2*quad_x-1,k);
                    
                    [avg1,jump1]=TraceVal1(2*[x1;xm;x2]-1,k);
                    
                    if l2==ratio*l1-1
                        jump1(3,:)=jump1(1,:);
                        jump1(1:2,:)=0;
                    end
                    if l2==min(ratio*(l1+1)-1,2^max(0,n2-1)-1)+1
                        jump1(1,:)=jump1(3,:);
                        jump1(2:3,:)=0;
                    end
                else
                    phi_1 = 2^((n1-1)/2)*waveletbasis(2^(n1-1)*quad_x-l1,k);
                    dphi_1 = 2^((n1-1)/2)*2^(n1-1)*dwaveletbasis(2^(n1-1)*quad_x-l1,k);
                    
                    [avg1,jump1]=TraceVal2(2^(n1-1)*[x1;xm;x2]-l1,k);
                    avg1=avg1*2^((n1)/2)*2^(n1-1);
                    jump1=jump1*2^((n1)/2);
                    
                    if l2==ratio*l1-1
                        % jump1
                        jump1(1:2,:)=0;
                        
                    end
                    if l2==min(ratio*(l1+1)-1,2^max(0,n2-1)-1)+1
                        
                        jump1(2:3,:)=0;
                        
                    end
                end
                
                
                if n2 == 0
                    if  l2>=ratio*l1 && l2<=min(ratio*(l1+1)-1,2^max(0,n2-1)-1)
                        phi_2 = legendre(2*quad_x-1,k);
                        dphi_2 = 2*dlegendre(2*quad_x-1,k);
                    else
                        phi_2 =zeros(2*quad_num,k);
                        dphi_2=zeros(2*quad_num,k);
                    end
                    
                    [avg2,jump2]=TraceVal1(2*[x1;xm;x2]-1,k);
                    if l2==ratio*l1-1
                        jump2(1:2,:)=0;
                    end
                    if l2==min(ratio*(l1+1)-1,2^max(0,n2-1)-1)+1
                        jump2(2:3,:)=0;
                    end
                else
                    if    l2>=ratio*l1 && l2<=min(ratio*(l1+1)-1,2^max(0,n2-1)-1)
                        phi_2=2^((n2-1)/2)*waveletbasis(2^(n2-1)*quad_x-l2,k);
                        dphi_2=2^((n2-1)/2)*2^(n2-1)*dwaveletbasis(2^(n2-1)*quad_x-l2,k);
                    else
                        phi_2 =zeros(2*quad_num,k);
                        dphi_2=zeros(2*quad_num,k);
                    end
                    
                    [avg2,jump2]=TraceVal2(2^(n2-1)*[x1;xm;x2]-l2,k);
                    avg2=avg2*2^((n2)/2)*2^(n2-1);
                    jump2=jump2*2^((n2)/2);
                    
                    if l2==ratio*l1-1
                        jump2(1:2,:)=0;
                    end
                    if l2==min(ratio*(l1+1)-1,2^max(0,n2-1)-1)+1
                        jump2(2:3,:)=0;
                    end
                end
                
                if l2==0
                    avg1(1,:)=2*avg1(1,:);
                    avg2(1,:)=2*avg2(1,:);
                end
                if l2==2^max(0,n2-1)-1
                    avg1(3,:)=2*avg1(3,:);
                    avg2(3,:)=2*avg2(3,:);
                end
                
                
                % (du,dv)
                dval=dphi_2'*(quad_w.*dphi_1);
                % <{u'},[v]>
                aujv=avg2'*jump1;
                % <{v'},[u]>
                avju=jump2'*avg1;
                %h^{-1}<[u],[v]>
                jujv=jump2'*jump1;
                
                % Generate the entry
                %                 val=dval-aujv-avju+sigma*jujv;
                if n1==n2 && l2==ratio*l1-1
                    jujv=zeros(k);
                end
                

                val=theta*dval-s1*aujv-s2*avju+sigma*jujv;
                
                mval=phi_1'*(quad_w.*phi_2);
                
                % Index of U and V
                if n1==0
                    Iu=[1:k];
                else
                    Iu=[k*2^(n1-1)+1+k*(l1):k*2^(n1-1)+k*(l1+1)];
                end
                if n2==0
                    Iv=[1:k];
                else
                    Iv=[k*2^(n2-1)+1+k*(l2):k*2^(n2-1)+k*(l2+1)];
                end

                [Iu,Iv]=meshgrid(Iu,Iv);
                
                if n1==n2 && l2<l1
                        val=zeros(k);     
                end
                if n1==n2 && l1==l2
                    M_stiff_diag=M_stiff_diag+sparse(Iu,Iv,val,dof,dof);
                    M_mass_diag=M_mass_diag+sparse(Iu,Iv,mval,dof,dof);
                else
                    
                    M_stiff=M_stiff+sparse(Iu,Iv,val,dof,dof);
                    M_mass=M_mass+sparse(Iu,Iv,mval,dof,dof);
                end

                
            end % loop of l2
        end % loop of n2
        
        
        
        % RHS
        x1=2^(-n1+1)*l1;
        if n1==0
            x2=1;
        else
            x2=2^(-n1+1)*(l1+1);
        end
        xm=(x1+x2)/2;

        % divide the whole to two different sub-cells
        
        [quad_x,quad_w]=lgwt(quad_num,x1,(x1+x2)/2);
        [quad_x2,quad_w2]=lgwt(quad_num,(x1+x2)/2,x2);
        
        quad_x=[quad_x;quad_x2];
        quad_w=[quad_w;quad_w2];
        

        if n1 ==0
            % generate the basis phi_1
            phi_1   = legendre(2*quad_x-1,k);
            E1 = legendre(2*xnode-1,k);
            
            Iu=[1:k];
            
        else
            phi_1 = 2^((n1-1)/2)*waveletbasis(2^(n1-1)*quad_x-l1,k);

            E1=2^((n1-1)/2)*waveletbasis(2^(n1-1)*xnode-l1,k);
            
            Iu=[k*2^(n1-1)+k*(l1)+1:k*2^(n1-1)+k*(l1+1)];
            
        end
        ff=rhs(quad_x);

        b(Iu)=phi_1'*(quad_w.*ff);
        
        tmp=exactu(quad_x);
        
        uu(Iu,1)=phi_1'*(quad_w.*tmp);

        % evalue basis
        
        Meval(:,Iu)=E1;
        
    end % loop of l1
end % loop of n1


M_stiff=M_stiff+M_stiff'+M_stiff_diag;


