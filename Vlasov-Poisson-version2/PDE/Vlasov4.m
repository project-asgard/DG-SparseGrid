function pde=Vlasov4
% Numerical Example for Vlasov Equation
% Bump-on-tail instability

% Parameters

k_0=0.3;
A=0.04;
Lmax=20*pi/3;
Vmax=13;

    % Initial condition 1D for x coordinate
    function f=Fx_0(x)
        
        f=(1+A*cos(k_0*x));
    end

    % Initial condition 1D for v coordinate
    function f=Fv_0(v)
        
        np=9/(10*sqrt(2*pi));
        nb=2/(10*sqrt(2*pi));
        u=4.5;
        vt=0.5;
        f=np*exp(-v.^2/2)+nb*exp(-(v-u).^2/(2*vt^2));
        
    end

    % Initial condition 2D
    function f=Fxv_0(x,v)
        
        f=Fv_0(v).*Fx_0(x);
        
    end

    % Electric field spatial variation
    function f=Ex(x)
        
        f=x-x;
        
    end

    % Electric field temporal variation
    function f=Et(x)
        
        f=x-x;
        
    end

    % Electric field
    function f=E(x,t)
        f=x-x;
    end

    function f=F(x,v,t)
        f=Fv_0(v).*(A*cos(k_0*(x-v*t)));
    end

    function f=rho(x,t)
        f= -A*cos(k_0*x).*exp(-k_0*t);
    end

pde = struct('k_0',k_0,'Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,'Ex',@Ex,...
    'Et',@Et,'E',@E,'F',@F,'rho',@rho,'Vmax',Vmax,'Lmax',Lmax);

end