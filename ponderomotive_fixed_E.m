function pde = ponderomotive_fixed_E
% Apply fixed E(x,t) and solve Vlasov only.
% E(x,t) has a gradient to allow ponderomotive force investigation. 

% parameters

A=1;
Lmax=1;
Vmax=5;
k_0=2*pi/(Lmax/2);

Lmin=0;Lmax=4*pi;
Vmin=-2*pi;Vmax=2*pi;

params.TEND = 1;

params.k_0 = k_0;
params.A = A;

params.Lmin = Lmin;
params.Lmax = Lmax;
params.Vmin = Vmin;
params.Vmax = Vmax;

pde.Fx_0 = @Fx_0;
pde.Fv_0 = @Fv_0;
pde.Fxv_0 = @Fxv_0;
pde.Ex = @Ex;
pde.Et = @Et;
pde.E = @E;
pde.rho = @rho;
pde.params = params;

    function f=Fx_0(x, params)
        % Initial condition for x variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
        
        offset = Lmax/2;
        
        lambda = 2*pi/k_0;
        xMin = offset-labmda/4;
        xMax = offset+lambda/4;
 
        f = x*0;
        here = (x>=xMin && x<=xMax);
        f(here)=(1+A*cos(k_0*(x(here)-offset)));
  
    end

    function f=Fv_0(v, params)
        % Initial condition for v variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        u=0.99;
        vt=0.3;
        f=(exp(-(u+v).^2/(2*vt^2))+exp(-(v-u).^2/(2*vt^2)))/(2*vt*sqrt(2*pi));
    end
    function f=Fxv_0(x,v, params)

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f=Fv_0(v).*Fx_0(x);
    end
    function f=Ex(x, params)

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f=x-x;
    end
    function f=Et(x, params)

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f=x-x;
        
    end
    function f=E(x,t,  params)

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
%         f=x-x;
%         f=A/k_0*exp(-k_0^2*t^2/2).*sin(k_0*x)-x+8*pi;
            f=x-x;%A/k_0*sin(k_0*x).*exp(-k_0^2*t.^2/2)-x+8*pi;
    end
    function f=F(x,v,t, params)

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f=Fv_0(v).*(A*cos(k_0*(x-v*t)));
    end
    function f=rho(x,t, params)

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f= -A*cos(k_0*x).*exp(-k_0*t);
    end
    
