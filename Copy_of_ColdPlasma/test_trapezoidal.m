% test for Trapezoidal rule

% d/dt G = F(G) with F=-G and G= exp(-t) 

G0 = 0;
CFL = 40;
dx = 2^(-10);
dt = CFL*dx;
% func = @(x)(-x);
 
MaxT=ceil(1/dt);

for t = 1:MaxT
% G1 = (G0-dt/2*(G0))/(1+dt/2);

% G1 = (G0+dt/2*(2));%/(1+dt/2); % t
G1 = (G0+dt*(2*t-1)*dt); % t^2

time = dt*t;
G0 = G1;
% [G1 time^2]
% [G1 exp(-time)]
end

max(abs(G1-time^2))