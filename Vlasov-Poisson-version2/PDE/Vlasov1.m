function pde=Vlasov1
% Numerical Example for Vlasov Equation
% Landau damping

% parameters
k_0=0.5; % k_0=2*pi/Lmax
A=0.5;%k_0;
Vmax=5;Lmax=4*pi;

params.k_0 = k_0;
params.A = A;
params.Vmax = Vmax;
params.Lmax = Lmax;


pde = struct('k_0',k_0,'Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,'Ex',@Ex,...
    'Et',@Et,'E',@E,'F',@F,'rho',@rho,'Vmax',Vmax,'Lmax',Lmax,  ...
    'params', params);
end

    function f=Fx_0(x, params)
        % Initial condition for x variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

%         f=(1+0.3*cos(k_0*x))/sqrt(2*pi);
        f=(1+A*cos(k_0*x));
    end

    function f=Fv_0(v, params)
        % Initial condition for v variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=exp(-v.^2/2)/sqrt(2*pi);
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


        % f=-k_0*sin(k_0*x)*exp(-k_0*t^2/2);
        f=-k_0*sin(k_0*x)
    end
    function f=Et(t, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=exp(-k_0*t^2/2);
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

        f=Fv_0(v).*(1+A*cos(k_0*(x-v*t)));
    end
    function f=rho(x,t, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

       f=1+exp(-k_0^2*t.^2/2)*A.*cos(k_0*x);
    end
    
