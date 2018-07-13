Dof_Hash = size(InvHash,2);
Dofs1 = Dof_Hash*Deg^3;
[M,N]=matrix_plot(4,2,1,FMWT_COMP_x);
s=size(M,2);
CoMat=sparse(s,Dofs1);
CoMat1=sparse(s^2,Dofs1);
for i = 1:Dof_Hash
    ll = InvHash{i};

    % Lev, Cell, Index from Hash
    iInd1 = ll(7);iInd2 = ll(8);iInd3 = ll(9);
    val=sparse(s,Deg^3);
    val1=sparse(s^2,Deg^3);
    count=1;
    for j = Deg*(iInd1-1)+1:Deg*iInd1
        for k= Deg*(iInd2-1)+1:Deg*iInd2
            for l= Deg*(iInd3-1)+1:Deg*iInd3
    %val3=M(:,Deg*(iInd3-1)+1:Deg*iInd3);
    val(:,count)=M(:,j).*M(:,k).*M(8,l);
    val1(:,count)=kron(M(:,j),M(:,k)).*M(8,l);
    count=count+1;
            end
        end
    end
    
    CoMat(:,Deg^3*(i-1)+1:Deg^3*i)=val;
    CoMat1(:,Deg^3*(i-1)+1:Deg^3*i)=val1;
end
E1h=Eh(1:Dofs1,1);
z0=-cos(2*pi*N).*sin(2*pi*N).*sin(2*pi*N(8)).*cos(pi/2*0.01);
z1=CoMat*E1h;
z2=CoMat1*E1h;
z3=reshape(z2,s,s);
[x,y]=meshgrid(N,N);
z=-cos(2*pi*x).*sin(2*pi*y).*sin(2*pi*N(8)).*cos(pi/2*0.01);

%Lev=4, fix z=2.3679e-01
