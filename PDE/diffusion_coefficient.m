v_par = 1:0.05:5;
v_perp = 1:0.05:5;
[V_PAR, V_PERP] = meshgrid(v_par,v_perp);
D_perp = zeros(1,length(V_PAR));
D_perp(V_PAR>=1 & V_PAR<2) = 1;
D_perp(V_PAR>=2 & V_PAR<3) = 1/2;
D_perp(V_PAR>=3 & V_PAR<4) = 1/3;
D_perp(V_PAR>=4 & V_PAR<5) = 1/4;
D_perp(V_PAR>=5) = 1/5;
%moving to spherical coordinates
v_mag = sqrt(V_PAR.^2 + V_PERP.^2);
theta = asin(V_PERP./v_mag);
D_perp = reshape(D_perp, [length(v_par) length(v_perp)]);
D_uu = (sin(theta)).^2.*D_perp;%cql3d transformation
%plotting coefficients
%contour(V_PAR,V_PERP,D_perp);
contour(v_mag, theta,D_uu);
xlim([0 7]);
ylim([0 pi/2]);
%plot(D_uu); 
%xlim([0 pi]);