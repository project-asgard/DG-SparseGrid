function f = Exa0(x,t)
global C E R

% a = 0.25;%3.1250e-02;
% b = 0.5;%6*a;
% f = x-x;
% ix = find((x-a).*(x-b)<0);
% f(ix) = 1;
% 
% f = exp(-t)*sin(pi*x);

% Analytical Solution
% % h0 = 3; h1 = 0.5; h2 = 1; h3 = 0.7; h4 = 3; h6 = 3;
% % f = (h0+...
% %      h1*x*exp(-1*(1+1)*t)+...
% %      h2*(3*x.^2-1)/2*exp(-2*(2+1)*t)+...
% %      h3*(5*x.^3-3*x)/2*exp(-3*(3+1)*t)+...
% %      h4*(35*x.^4-30*x.^2+3)/8*exp(-4*(4+1)*t)+...
% %      h6*(231*x.^6-315*x.^4+105*x.^2-5)/16*exp(-6*(6+1)*t));

% E = 2; C = 1;
% A = E/C;
% % f = A/(2*sinh(A))*exp(A*x);
% % f = A-x+x;
% f = (exp(1)-exp(-1))/sinh(A);

% E = 2;
% C = 1;
% R = 2;
A = E/C;
B = R/C;
% % % Q = 1;
[quad_x,quad_w]=lgwt(100,-1,1);
Q = sum(quad_w.*exp(A*quad_x+(B/2)*quad_x.^2))/2;
f = Q-x+x;
f = exp(A*x+(B/2)*x.^2);

% f = x-x+1;

% Test 3_2_3b
% sigma = 0.1;
% % f = exp(-x.^2/sigma^2);
% 
% % 
% f = exp(-(x+0.36).^2/sigma^2)+exp(-(x-0.36).^2/sigma^2);

a = 1;
f = 4/sqrt(pi)/a^3*exp(-x.^2/a^2);
% [quad_x,quad_w]=lgwt(100,0,10);
% Q = sum(quad_w.*4/sqrt(pi).*exp(-quad_x.^2/a^2));
% f = 4/sqrt(pi)*exp(-x.^2/a^2)/Q;

f = x-x+3/5^3;
ix = find(x>5);
f(ix) = 0;

end