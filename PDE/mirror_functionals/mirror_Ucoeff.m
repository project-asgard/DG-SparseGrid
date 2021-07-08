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
N = 100;
z = linspace(0,pi,N);
% y = myfun(1,z);
% y2 = myfun(1,z');
% plot(z,y);
U_coeff = zeros(size(uVals));
for i = 1:numel(uVals)
    U_coeff(i) = trapz(z,(2*l+1)/2.*myfun(uVals(i),z,l));
    %U_coeff(i) = integral(@(x) (2*l + 1)/2.*myfun(uVals(i),x),0,pi);
end
end