function [RHS_MWDG]=RHS_1Dv2(n,k,t,pde)
% for every time step t, compute the rhs
%--DG parameters
quad_num=10;
%---------------

% compute the trace values
p_1 = legendre(-1,k);
p_2 = legendre(1,k);

[quad_x,quad_w]=lgwt(quad_num,-1,1);

p_val = legendre(quad_x,k);

dof_1D=k*2^n;
b3x1=sparse(dof_1D,1);
b3x2=sparse(dof_1D,1);
b3x3=sparse(dof_1D,1);
b3y1=sparse(dof_1D,1);
b3y2=sparse(dof_1D,1);
b3y3=sparse(dof_1D,1);
b3y1=sparse(dof_1D,1);
b3y2=sparse(dof_1D,1);
b3y3=sparse(dof_1D,1);

b4x1=sparse(dof_1D,1);
b4x2=sparse(dof_1D,1);
b4x3=sparse(dof_1D,1);
b4y1=sparse(dof_1D,1);
b4y2=sparse(dof_1D,1);
b4y3=sparse(dof_1D,1);
b4y1=sparse(dof_1D,1);
b4y2=sparse(dof_1D,1);
b4y3=sparse(dof_1D,1);

% generate 1D rhs for DG 
for LL=0:2^n-1
    % RHS term
     f3x=pde.rhs3x( quad_x*2^(-n-1)+2^(-n)*(LL+0.5) );
     f3y=pde.rhs3y( quad_x*2^(-n-1)+2^(-n)*(LL+0.5) );
     f3z=pde.rhs3z( quad_x*2^(-n-1)+2^(-n)*(LL+0.5) );
     
     f4x=pde.rhs4x( quad_x*2^(-n-1)+2^(-n)*(LL+0.5) );
     f4y=pde.rhs4y( quad_x*2^(-n-1)+2^(-n)*(LL+0.5) );
     f4z=pde.rhs4z( quad_x*2^(-n-1)+2^(-n)*(LL+0.5) );
     
     Iu=[k*LL+1:k*(LL+1)];
     b3x1(Iu)=p_val'*(quad_w.*f3x(:,1))*2^(-n)/2*2^(n/2);
     b3x2(Iu)=p_val'*(quad_w.*f3x(:,2))*2^(-n)/2*2^(n/2);
     b3x3(Iu)=p_val'*(quad_w.*f3x(:,3))*2^(-n)/2*2^(n/2);
     b3y1(Iu)=p_val'*(quad_w.*f3y(:,1))*2^(-n)/2*2^(n/2);
     b3y2(Iu)=p_val'*(quad_w.*f3y(:,2))*2^(-n)/2*2^(n/2);
     b3y3(Iu)=p_val'*(quad_w.*f3y(:,3))*2^(-n)/2*2^(n/2);
     b3z1(Iu)=p_val'*(quad_w.*f3z(:,1))*2^(-n)/2*2^(n/2);
     b3z2(Iu)=p_val'*(quad_w.*f3z(:,2))*2^(-n)/2*2^(n/2);
     b3z3(Iu)=p_val'*(quad_w.*f3z(:,3))*2^(-n)/2*2^(n/2);
     
     b4x1(Iu)=p_val'*(quad_w.*f4x(:,1))*2^(-n)/2*2^(n/2);
     b4x2(Iu)=p_val'*(quad_w.*f4x(:,2))*2^(-n)/2*2^(n/2);
     b4x3(Iu)=p_val'*(quad_w.*f4x(:,3))*2^(-n)/2*2^(n/2);
     b4y1(Iu)=p_val'*(quad_w.*f4y(:,1))*2^(-n)/2*2^(n/2);
     b4y2(Iu)=p_val'*(quad_w.*f4y(:,2))*2^(-n)/2*2^(n/2);
     b4y3(Iu)=p_val'*(quad_w.*f4y(:,3))*2^(-n)/2*2^(n/2);
     b4z1(Iu)=p_val'*(quad_w.*f4z(:,1))*2^(-n)/2*2^(n/2);
     b4z2(Iu)=p_val'*(quad_w.*f4z(:,2))*2^(-n)/2*2^(n/2);
     b4z3(Iu)=p_val'*(quad_w.*f4z(:,3))*2^(-n)/2*2^(n/2);
     
%      b(Iu)=p_val'*(quad_w.*f3z(:,1))*2^(-n)/2*2^(n/2);
     
%      
%      
%      Meval(2*LL+1:2*LL+2,Iu)=2^(n/2)*legendre([-1,1],k);

end
