function U_coeff = mirror_Ucoeff(f_prev,uVals,l)
    function ans = myfun(u,z)
        assert(numel(u)==1);
        nz = numel(z);
        sizez = size(z);
        z = reshape(z,1,nz);
        legendre_vals = legendre(l, z);
        legendre_val = legendre_vals(l+1,:);
        f = f_prev(u,z);
        ans = f.*legendre_val;
        ans = reshape(ans,sizez);
    end
% z = linspace(-1,1);
% y = myfun(1,z);
% y2 = myfun(1,z');
% plot(z,y);
for i = 1:numel(uVals)
    U_coeff(i) = integral(@(z) myfun(uVals(i),z), -1, 1);
end
end