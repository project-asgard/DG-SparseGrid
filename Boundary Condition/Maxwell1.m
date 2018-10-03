function pde=Maxwell1
%====================================
% Numerical Test 1 for Maxwells' Eq.
% Non-periodic boundary condition
%====================================


c0=3e8;
w=pi/2;%2*pi*sqrt(3)/2*c0;

    function E=E(x)
        E= cos(pi*x);
    end

    function f=f(x)
        f= -pi*sin(pi*x);
    end

pde = struct('w',w,'E',@E,'f',@f);

end