function pde=Vlasov5
% Numerical Example for Vlasov Equation
% Two-stream instability

% parameters
k_0=0.5;A=0.05;Lmax=4*pi;Vmax=2*pi;

Lmin=0;Lmax=4*pi;
Vmin=-2*pi;Vmax=2*pi;


% params.k_0 = k_0;
% params.A = A;
% params.Lmax = Lmax;
% params.Vmax = Vmax;

% pde.Fx_0 = @Fx_0;
% pde.Fv_0 = @Fv_0;
% pde.Fxv_0 = @Fxv_0;
% pde.Ex = @Ex;
% pde.Et = @Et;
% pde.E = @E;
% pde.rho = @rho;
% pde.params = params;

% Modified By Lin
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

pde.solvePoisson = 1;
pde.applySpecifiedE = 0;
pde.implicit = 0;

pde.checkAnalytic = 0;

pde.exactE = @exactE;

pde.source1x = @source1x;
pde.source1v = @source1v;
pde.source1t = @source1t;
pde.source1 = @source1;

pde.source2x = @source2x;
pde.source2v = @source2v;
pde.source2t = @source2t;
pde.source2 = @source2;

pde.source3x = @source3x;
pde.source3v = @source3v;
pde.source3t = @source3t;
pde.source3 = @source3;

pde.ExactFx = @ExactFx;
pde.ExactFv = @ExactFv;
pde.ExactFt = @ExactFt;

% pde = struct('k_0',k_0,'Fx_0',@Fx_0,'Fv_0',@Fv_0,'Fxv_0',@Fxv_0,'Ex',@Ex,...
%     'Et',@Et,'E',@E,'F',@F,'rho',@rho,'Vmax',Vmax,'Lmax',Lmax);

end


    function f=Fx_0(x, params)
        % Initial condition for x variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f=(1+A*cos(k_0*x));
    end

    function f=Fv_0(v, params)
        % Initial condition for v variable

        A = params.A;
        k_0 = params.k_0;
        Lmax = params.Lmax;
        Vmax = params.Vmax;
 
        f=v.^2.*exp(-v.^2/2)/sqrt(2*pi);
    end
    function f=Fxv_0(x,v,  params)

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
 
        f=x-x;
    end
    function f=Et(x,  params)

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
    function f=F(x,v,t,  params)

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
 
%        f=exp(-k_0*t.^2/2).*cos(k_0*x); 
        f= -A*cos(k_0*x).*exp(-k_0*t);
    end
    
    % Added by Lin
    function f = exactE(x, params)
% Exact solution for E
f=x*0;
end

% Source term 1
function f = source1x(x, params)
f = x*0;
end

function f = source1t(t)
f = t*0;
end

function f = source1v(v, params)
f = v*0;
end

function f = source1(x,v,t)
f = source1x(x).*source1v(v).*source1t(t);
end

% Source term 2
function f = source2x(x, params)
f = x*0;
end

function f = source2v(v, params)
f = v*0;
end

function f = source2t(t)
f = t*0;
end

function f = source2(x,v,t)
f = source2x(x).*source2v(v).*source2t(t);
end

% source term 3
function f = source3x(x, params)
f = x*0;
end

function f = source3v(v, params)
f = v*0;
end

function f = source3t(t)
f = t*0;
end

function f = source3(x,v,t)
f = source3x(x).*source3v(v).*source3t(t);
end

% Exact F
function f=ExactFx(x, params)
f = x*0;
end

function f=ExactFv(v, params)
f = v*0;
end

function f=ExactFt(t)
f = t*0;
end

function f=ExactF(x,v,t)
f = ExactFx(x).*ExactFv(v).*ExactFt(t);
end
