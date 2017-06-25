% 1D approximation for Poisson Eq

dof = 2^n*k;

[Stiff_1D,b,xnode,Meval,uu] = LaplacianMatrix(n,k);
M_mass = diag(ones(1,dof));
sol_1D = Stiff_1D\b*pi^2;

%--------------
% 1D plot
%--------------
figure
plot(xnode,Meval*sol_1D,'r+');
hold on
plot(xnode,exactu(xnode),'b')
plot(xnode,Meval*uu,'go')
legend('Numerical Sol','Tru Sol','u_I')


Max1D=max(abs(uu(2:end)-sol_1D(2:end)));
disp(['Max Error of 1D Computing   ' num2str(Max1D)])