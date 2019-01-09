% Test for analytic solution
E = 1;
C = 1;
R = 1;

ExactF = @(x,t)(sin(pi*x)*t);
source = @(x,t)( sin(pi*x)+...
    +E*t*(-2*x.*sin(pi*x)+(-x.^2+1).*cos(pi*x)*pi)+...
    -C*t*(-2*x.*cos(pi*x)*pi-(-x.^2+1).*sin(pi*x)*pi^2)+...
    +R*t*( (-x.^2+1).*sin(pi*x)-2*x.^2.*sin(pi*x)+x.*(-x.^2+1).*cos(pi*x)*pi) );
