function pde=Maxwell1
%====================================
% Numerical Test 1 for Maxwells' Eq.
%====================================


c0=3e8;
w=pi/2;%2*pi*sqrt(3)/2*c0;

    function E=exact_E1(x,t)
        E= -cos(2*pi*x).*cos(w*t);
    end

    function E=E1(x)% without t
        E= -cos(2*pi*x);
    end
    
    function E=exact_E2(x,t)
        E= sin(2*pi*x).*cos(w*t);
    end

    function E=E2(x)
        E= sin(2*pi*x);
    end
    
    function B=exact_B(x,t)
        
        B= 2*pi/w*sin(2*pi*x).*sin(w*t);
    end

    function B=B(x)
        B= 2*pi/w*sin(2*pi*x);
    end

    function f=exact_f1(x,t)
        f= (4*pi*pi/w-w)*cos(2*pi*x).*sin(w*t);
    end

    function f=f1(x)
        f= (4*pi*pi/w-w)*cos(2*pi*x);
    end

    function f=exact_f2(x,t)
        f= -w*sin(2*pi*x).*sin(w*t);
    end

    function f=f2(x)
        f= -w*sin(2*pi*x);
    end

pde = struct('w',w,'E1',@exact_E1,'E2',@exact_E2,'B',@exact_B,'F1',@exact_f1,'F2',@exact_f2,...
    'f1',@f1,'f2',@f2,'e1',@E1,'e2',@E2,'b',@B);

end