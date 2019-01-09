%This code presents the analytical solution of the collision
h0 = 3; h1 = 0.5; h2 = 1; h3 = 0.7; h4 = 3; h6 = 3;
ExactF = @(x,t)(h0+h1*x+h2*(3*x.^2-1)/2+h3*(5*x.^3-3*x)/2+...
    h4*(35*x.^4-30*x.^2+3)/8+h6*(231*x.^6-315*x.^4+105*x.^2-5)/16);
source = @(x,t)(x-x);


E = 0;
C = 1;
R = 0;