function f = TimeAdvance(A,f,dt)
%-------------------------------------------------
% Time Advance Method
% Input: Matrix:: A
%        Vector:: f
%        Time Step:: dt
% Output: Vector:: f
%-------------------------------------------------

f = RungeKutta3(A,f,dt);

end

function fval = RungeKutta3(A,f,dt)
%----------------------------------
% 3-rd Order Runge Kutta Method
%----------------------------------
ftmp = ApplyA(A,f);
f_1 = f +dt *ftmp;
ftmp = ApplyA(A,f_1);
f_2 = 3/4*f+1/4*f_1+1/4*dt*(ftmp);
ftmp = ApplyA(A,f_2);
fval = 1/3*f+2/3*f_2+2/3*dt*(ftmp);

end

function ftmp = ApplyA(A,f)
%-----------------------------------
% Multiply Matrix A by Vector f
%-----------------------------------
dof = size(f,1);
ftmp=sparse(dof,1);
for i=1:size(A,2)
tmpA=A{i}.A1;
tmpB=A{i}.A2;
IndexI=A{i}.IndexI;
IndexJ=A{i}.IndexJ;
ftmp(IndexI)=ftmp(IndexI)+kron_mult2(tmpA,tmpB,f(IndexJ));
end

end







