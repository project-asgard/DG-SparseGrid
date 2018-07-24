% plot for numerical and exact solutions
M_2Deval = kron(Meval,Meval);

[x_2D,y_2D] = meshgrid(xnode,xnode);

figure(10);

subplot(1,3,1)
mesh(x_2D,y_2D,reshape(M_2Deval*sol_2D,size(x_2D)));
title('Full Grid Solution')
axis([0 1 0 1 0 1 ])


subplot(1,3,2)
mesh(x_2D,y_2D,reshape(M_2Deval*Iu_2D,size(x_2D)));
title('Exact Solution')
axis([0 1 0 1 0 1 ])


subplot(1,3,3)
mesh(x_2D,y_2D,reshape(M_2Deval*(S2F'*sol_s),size(x_2D)));
title('Sparse Grid Approximation')
axis([0 1 0 1 0 1 ])