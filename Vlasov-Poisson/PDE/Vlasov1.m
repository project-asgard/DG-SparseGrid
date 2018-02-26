function pde=Vlasov1
% Numerical Example for Vlasov Equation
% Landau damping

% parameters
k_0=0.5; % k_0=2*pi/Lmax
A=0.5;%k_0;
Vmax=5;Lmax=4*pi;


    function f=Fx_0(x)
        % Initial condition for x variable
%         f=(1+0.3*cos(k_0*x))/sqrt(2*pi);
        f=(1+A*cos(k_0*x));
    end

    function f=Fv_0(v)
        % Initial condition for v variable
        f=exp(-v.^2/2)/sqrt(2*pi);
    end
    function f=Fxv_0(x,v)
        f=Fv_0(v).*Fx_0(x);
    end
    function f=Ex(x)
        f=-k_0*sin(k_0*x)*exp(-k_0*t^2/2);
    end
    function f=Et(t)
        f=exp(-k_0*t^2/2);
    end
    function f=E(x)
        f=Ex(x).*Et(t);
    end
    function f=F(x,v,t)
        f=Fv_0(v).*(1+A*cos(k_0*(x-v*t)));
    end
    function f=rho(x,t)
       f=1+exp(-k_0^2*t.^2/2)*A.*cos(k_0*x);
    end
    
pde = struct('k_0',k_0,'Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,'Ex',@Ex,...
    'Et',@Et,'E',@E,'F',@F,'rho',@rho,'Vmax',Vmax,'Lmax',Lmax);
end