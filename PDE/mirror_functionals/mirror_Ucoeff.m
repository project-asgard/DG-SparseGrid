function V_coeff = mirror_Ucoeff(f_prev,lIndex,uval)

    dz = 5e-5;
    num_cells = pi/dz;
    V_coeff = 0;
    for k = 0:num_cells
        z = -pi + k*dz;
        val_z = cos(z);
        val_z1 = cos(z+dz);
        legendre_Val = legendre(lIndex,val_z);
        legendre_Val1 = legendre(lIndex,val_z1);
        integrand = ((-f_prev(uval,z)*legendre_Val(lIndex + 1)*sin(z))+ (-f_prev(uval,z+dz)*legendre_Val1(lIndex + 1)*sin(z+dz)))*dz/2;
        V_coeff = V_coeff + integrand;
    end
    V_coeff = (2*lIndex + 1)*V_coeff/2;
end