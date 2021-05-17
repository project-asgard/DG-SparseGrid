function func = mirror_background(f_prev,lmax,uval)

    func = @(u,z) u.*0;
    for l = 1:lmax
        integrand = f_prev(uval,z).*legendreP(l,z);
        func(uval,z) = func(uval,z) + mirror_V_coeff(f_prev,l,uval).*legendreP(l,z);
    end

end
