function pde=Vlasov4
% Numerical Example for Vlasov Equation
% Bump-on-tail instability

% Parameters

k_0=0.3;
A=0.04;
Lmax=20*pi/3;
Vmax=13;

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
pde.params = params;

% pde = struct('k_0',k_0,'Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,'Ex',@Ex,...
%     'Et',@Et,'E',@E,'F',@F,'rho',@rho,'Vmax',Vmax,'Lmax',Lmax, ...
%     'params',    params);

end

    % Initial condition 1D for x coordinate
    function f=Fx_0(x,  params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f=(1+A*cos(k_0*x));
    end

    % Initial condition 1D for v coordinate
    function f=Fv_0(v,  params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
         
        np=9/(10*sqrt(2*pi));
        nb=2/(10*sqrt(2*pi));
        u=4.5;
        vt=0.5;
        f=np*exp(-v.^2/2)+nb*exp(-(v-u).^2/(2*vt^2));
        
    end

    % Initial condition 2D
    function f=Fxv_0(x,v, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
        
        f=Fv_0(v).*Fx_0(x);
        
    end

    % Electric field spatial variation
    function f=Ex(x, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
        
        f=x-x;
        
    end

    % Electric field temporal variation
    function f=Et(x, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
        
        f=x-x;
        
    end

    % Electric field
    function f=E(x,t,  params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f=x-x;
    end

    function f=F(x,v,t, params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;


        f=Fv_0(v).*(A*cos(k_0*(x-v*t)));
    end

    function f=rho(x,t,  params)
        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;

        f= -A*cos(k_0*x).*exp(-k_0*t);
    end

