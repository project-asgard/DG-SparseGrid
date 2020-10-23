v_par = -5:0.05:5;
v_per = 0:0.05:5;
[V_PAR, V_PER] = meshgrid(v_par,v_per);
D_per = V_PAR*0;
step_indices = V_PAR>1 & V_PAR<3;
D_per(step_indices(:))=1;

%moving to spherical coordinates
v_mag = sqrt(V_PAR.^2 + V_PER.^2);
theta = asin(V_PER./v_mag);
D_uu = (sin(theta)).^2.*D_per;%cql3d transformation

%plotting coefficients
subplot(2,1,1)
contourf(V_PAR,V_PER,D_per);
subplot(2,1,2)
contourf(v_mag,theta,D_uu);
xlim([0 7]);
ylim([0 pi/2]);
%plot(D_uu); 
%xlim([0 pi]);