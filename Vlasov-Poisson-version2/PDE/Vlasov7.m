function pde=Vlasov7
% Numerical Example for Vlasov Equation
% This test has given E
% need to check whether E is correct

% parameters
k_0=0.5; 
A=1;
Vmax=1;Lmax=1;

params.k_0 = k_0;
params.A = A;
params.Lmax = Lmax;
params.Vmax = Vmax;

pde.Fx_0 = @Fx_0;
pde.Fv_0 = @Fv_0;
pde.Fxv_0 = @Fxv_0;
pde.Ex = @Ex;
pde.Et = @Et;
pde.E = @E;
pde.rho = @rho;
pde.IsExactE = IsExactE;
pde.exactE = @exactE;
pde.params = params;

% pde = struct('k_0',k_0,'Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,'Ex',@Ex,...
%     'Et',@Et,'E',@E,'F',@F,'rho',@rho,'Vmax',Vmax,'Lmax',Lmax);
end

    function f=Fx_0(x, params)
        % Initial condition for x variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=cos(2*pi*x);
    end

    function f=Fv_0(v, params)
        % Initial condition for v variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=exp(-k_0*v.^2/2);
    end
    function f=Fxv_0(x,v, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=Fv_0(v).*Fx_0(x);
    end
    function f=Ex(x,  params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=-2*pi/k_0*tan(2*pi*x);
    end
    function f=Et(t, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=t-t+1;
    end
    function f=E(x,t, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=Ex(x).*Et(t);
    end
    function f=F(x,v,t,  params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=exp(-t)*exp(-k_0*v.^2/2).*cos(2*pi*x);
    end
    function f=rho(x,t, params)
    % we do note use rho in this test, so I just write a arbitrary function
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

       f=x-x+1;
    end
    function f=IsExactE
        % IsExactE = 1 means given function for E

        f=1;
    end
    function f = exactE(x)
    % Exact solution for E
        k_0 = 0.5;
        f=-2*pi/k_0*tan(2*pi*x);
    end
    function f = sourse(x,v,t)
    % source term
        f = 
    end

    
