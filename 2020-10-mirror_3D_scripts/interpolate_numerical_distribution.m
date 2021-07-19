function F = interpolate_numerical_distribution()

pitch_file = fopen('pitch_i.txt', 'r');
vel_file = fopen('u_j.txt', 'r');
dist_file = fopen('f_ij.txt', 'r');

pitch_vals = fscanf(pitch_file, '%f');
vel_vals = fscanf(vel_file, '%f');
f_vals = fscanf(dist_file, '%f');

size_f = length(f_vals);
num = sqrt(size_f);

n_o = 8e14;
temp = 1.6022e-10; %temperature in erg
m_e = 9.109*10^-28;
u_th = sqrt(2*temp/m_e);
lIndex = [0 1 2];

[u,z] = meshgrid(vel_vals,pitch_vals);
maxwell = n_o.*exp(-u.^2./u_th^2).*(z.*0 + 1)./(u_th.^3.*pi^(3/2));

f_vals = reshape(f_vals, [num num]);

ax1 = subplot(2,2,1);
mesh(vel_vals,pitch_vals,f_vals);
ax2 = subplot(2,2,2);
mesh(u,z,maxwell);

u = u';

z = z';

F = griddedInterpolant(u,z,f_vals);


end