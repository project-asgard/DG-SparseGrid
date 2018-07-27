function f = ImplicitTime2(s1,s2,h1,h2,A,f,dt,b1,b2) %Gauss-Legendre
%Still have some problem
%-------------------------------------------------
% Time Advance Method
% Input: Matrix:: A
%        Vector:: f
%        Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------
f = Implicit(s1,s2,h1,h2,A,f,dt,b1,b2);

end

function fval = Implicit(s1,s2,h1,h2,A,f,dt,b1,b2)
%----------------------------------
% Implicit Method
%----------------------------------
tmp1=s1*(b2-b1);
tmp2=s2*(b1-b2);
b1=(1/4-sqrt(3)/6)*A*dt*tmp1+b1;
b2=(1/4+sqrt(3)/6)*A*dt*tmp2+b2;
k1=h1*(A*f+b1);
k2=h2*(A*f+b2);
fval=f+1/2*k1*dt+1/2*k2*dt;
%ftmp=H*(A*f+b);
%fval=f+dt*ftmp;
%     sol_1=sol_n+dt*(A_s*sol_n+b_s*sin(pde.w*time));
%     sol_2=3/4*sol_n+1/4*sol_1+1/4*dt*(A_s*sol_1+b_s*sin(pde.w*time));
%     sol_n=1/3*sol_n+2/3*sol_2+2/3*dt*(A_s*sol_2+b_s*sin(pde.w*time));
end

function ftmp = ApplyA(A,f)
%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
ftmp=sparse(dof,1);

ftmp = A*f;
% for i=1:size(A,2)
% tmpA=A{i}.A1;
% tmpB=A{i}.A2;
% IndexI=A{i}.IndexI;
% IndexJ=A{i}.IndexJ;
% ftmp(IndexI)=ftmp(IndexI)+kron_mult2(tmpA,tmpB,f(IndexJ));
% end

end







