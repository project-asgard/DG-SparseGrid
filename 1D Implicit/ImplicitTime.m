function f = ImplicitTime(A,f,dt,b) %Backward Euler
%-------------------------------------------------
% Time Advance Method
% Input: Matrix:: A
%        Vector:: f
%        Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------
f = Implicit(A,f,dt,b);

end

function fval = Implicit(A,f,dt,b)
%----------------------------------
% Implicit Method
%----------------------------------

%condest(S);
%fval=H*(1/dt)*f+H*b;
fval=A*(f+b*dt);

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







