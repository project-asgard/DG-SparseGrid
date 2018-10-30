function pde=Maxwell13
%====================================
% Numerical Test 1 for Maxwells' Eq.
% Non-periodic boundary condition
%====================================


c0=3e8;
w=pi/2;%2*pi*sqrt(3)/2*c0;

    function E=exact_E1(x,t)
        E= -cos(pi*x).*cos(w*t);
    end

    function E=E1(x)% without t
        E= -cos(pi*x);
    end
    
    function E=exact_E2(x,t)
        E= sin(pi*x).*cos(w*t);
    end

    function E=E2(x)
        E= sin(pi*x);
    end
    
    function B=exact_B(x,t)
        
        B= pi/w*sin(pi*x).*sin(w*t);
    end

    function B=B(x)
        B= pi/w*sin(pi*x);
    end

    function f=exact_f1(x,t)
        f= (pi*pi/w-w)*cos(pi*x).*sin(w*t);
    end

    function f=f1(x)
        f= (pi*pi/w-w)*cos(pi*x);
    end

    function f=exact_f2(x,t)
        f= w*sin(pi*x).*sin(w*t);
    end

    function f=f2(x)
        f= w*sin(pi*x);
    end

    function f=exact_f3(x,t)
        f= x-x;
    end

    function f=f3(x)
        f= x-x;
    end

pde = struct('w',w,'E1',@exact_E1,'E2',@exact_E2,'B',@exact_B,'F1',@exact_f1,'F2',@exact_f2,...
    'f1',@f1,'f2',@f2,'e1',@E1,'e2',@E2,'b',@B,'F3',@exact_f3,'f3',@f3);

end