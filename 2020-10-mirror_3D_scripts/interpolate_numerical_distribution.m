function F_out = interpolate_numerical_distribution()

    function out = F2(u,z)
        out = ones(numel(u),numel(z));
        if size(u) == size(z)
            for i = 1:numel(u)
                for j = 1:numel(z)
                    out(i,j) = F(u(i),z(j));
                end
            end
        elseif isscalar(u) 
            for j = 1:numel(z)
                out(j) = F(u,z(j));
            end
        end
    end

pitch_file = fopen('pitch_i.txt', 'r');
vel_file = fopen('u_j.txt', 'r');
dist_file = fopen('f_ij.txt', 'r');

z = fscanf(pitch_file, '%f');
u = fscanf(vel_file, '%f');
f_vals = fscanf(dist_file, '%f');

size_f = length(f_vals);
num = sqrt(size_f);

n_o = 8e14;
temp = 1.6022e-10; %temperature in erg
m_e = 9.109*10^-28;
u_th = sqrt(2*temp/m_e);
lIndex = [0 1 2];

%[u,z] = meshgrid(vel_vals,pitch_vals);
%maxwell = n_o.*exp(-u.^2./u_th^2).*(z.*0 + 1)./(u_th.^3.*pi^(3/2));

analytic = @(x,y) n_o.*exp(-x.^2./u_th^2).*(y.*0 + 1)./(u_th.^3.*pi^(3/2));
f_vals = reshape(f_vals, [num num]);

% ax1 = subplot(1,3,1);
% contour(vel_vals,pitch_vals,f_vals);
% ax2 = subplot(1,3,2);
% contour(u,z,analytic(u,z));

%u = u';

%z = z';

F = griddedInterpolant({u,z},f_vals');
F_out = @F2;

%ax3 = subplot(1,3,3);
%f_plot = F_out(u(:),z(:));
%f_plot = reshape(f_plot, size(u));
%contour(u,z,f_plot);
% 



end