%This code presents the solution for Full Pitch Angle

E = 1; C = 1; R = 3;
A = E/C;
B = R/C;
Q = 1;
ExactF = @(x,t)(Q*exp(A*x+B/2*x.^2));
source = @(x,t)(x-x);