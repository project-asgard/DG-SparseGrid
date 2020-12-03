function U_coeff = mirror_Ucoeff(f_prev,uVals,l)
global abs_tol
global rel_tol
    function ans = myfun(u,z,l)
        assert(numel(u)==1);
        nz = numel(z);
        sizez = size(z);
        z = reshape(z,1,nz);
        legendre_vals = legendre(l, cos(z));
        legendre_val = legendre_vals(1,:);
        f = f_prev(u,z).*sin(z);
        ans = f.*legendre_val;
        ans = reshape(ans,sizez);
    end
pitch_file = fopen('pitch_long_i.txt', 'r');
dpitch_file = fopen('dpitch_long_i.txt', 'r');

%z = fscanf(pitch_file, '%f');
%dz = fscanf(dpitch_file, '%f');
N = 300;

% for i = 1:numel(z)-1
%     z(i) = (z(i+1) + z(i))/2;
% end
% DEFINE: Variable Density Vector — 
% s: 3-element vector of start values
% e: 3-element vector of end values
% d: 3-element vector of density values (number of entries between start and end) 
vdv = @(s,e,d) sort([linspace(s(1),e(1),d(1))  linspace(s(2),e(2),d(2))  linspace(s(3),e(3),d(3))]);
s = [0 pi/3+0.1 2*pi/3+0.1];
e = [pi/3 2*pi/3 pi];
d = [100 100 100];
%z = vdv(s,e,d);
% for i = 1:numel(z)-1
%     dz(i+1) = z(i+1) - z(i);
% end
% dz(1) = 0;
z = linspace(0,pi,N);
dz = pi/N;
% y = myfun(1,z);
% y2 = myfun(1,z');
% plot(z,y);
U_coeff = zeros(size(uVals));
for i = 1:numel(uVals)
    U_coeff(i) = (2.*l+1)./2.*sum(myfun(uVals(i),z,l).*dz);
    %U_coeff(i) = (2.*l+1)./2.*trapz(z,myfun(uVals(i),z,l));
    %U_coeff(i) = integral(@(x) (2*l + 1)/2.*myfun(uVals(i),x,l),0,pi);
end
end