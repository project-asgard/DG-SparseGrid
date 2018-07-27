function pde=Maxwell10
%====================================
% Numerical Test 1 for Maxwells' Eq.
%====================================


c0=3e8;
w=pi/2;%2*pi*sqrt(3)/2*c0;

    function E=exact_E1(x,t)
        E= t;
    end

    function E=E1(x)% without t
        E= 1;
    end
    
    function E=exact_E2(x,t)
        E= t;
    end

    function E=E2(x)
        E= 1;
    end
    
    function B=exact_B(x,t)
        
        B= t;
    end

    function B=B(x)
        B= 1;
    end

    function f=exact_f1(x,t)
        f= -1;
    end

    function f=f1(x)
        f= -1;
    end

    function f=exact_f2(x,t)
        f= -1;
    end

    function f=f2(x)
        f= -1;
    end

    function f=exact_f3(x,t)
        f= -1;
    end

    function f=f3(x)
        f= -1;
    end

pde = struct('w',w,'E1',@exact_E1,'E2',@exact_E2,'B',@exact_B,'F1',@exact_f1,'F2',@exact_f2,...
    'f1',@f1,'f2',@f2,'e1',@E1,'e2',@E2,'b',@B,'F3',@exact_f3,'f3',@f3);

end