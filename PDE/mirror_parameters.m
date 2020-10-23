function p = mirror_parameters()

% Common constants and functions for the mirror set of PDEs

m_e = 9.109*10^-31; %electron mass in kg
m_D = 3.3443*10^-27; %deuterium mass in kgs
m_H = 1.6726*10^-27; %hydrogen mass in kgs
k_b = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
e = 1.602*10^-19; %charge in Coulombs
ln_delt = 10; %Coulomb logarithm
v_th = @(T_eV,m) sqrt(2*T_eV * e/m);
eps0 = 8.85*10^-12; %permittivity of free space in Farad/m
mu0 = 4*pi*10^-7; %magnetic permeability in Henries/m
n_turns = 5; %number of turns in the mirror
R_mag = 0.4; %radius of individual loop in meters
I_mag = 10; %current going through individual loop in meters
B_o = mu0*n_turns*I_mag/(2*R_mag); %magnetic field under current loop
v_test = 500; %test value for velocity in coefficient
z_test = pi/2 - 1e-6; %test value for pitch angle in coefficient
advec_space_1D = v_test*cos(z_test); %advection coefficient for 1D in space

% species b: electrons in background
b.n = 4e19;
b.T_eV = 4;
b.Z = 1;
b.m = m_e;
b.vth = v_th(b.T_eV,b.m);

%species b2: deuterium in background
b2.n = 4e19;
b2.T_eV = 4;
b2.Z = 1;
b2.m = m_D;
b2.vth = v_th(b2.T_eV,b2.m);

% species a: electrons in beam
a.n = 4e19;
a.T_eV = 1e3;
a.Z = -1;
a.m = m_e;
a.vth = v_th(a.T_eV,a.m);

x       = @(v,vth) v./vth; %normalized velocity to thermal velocity
xi      = @(s,R_mag) s/R_mag; %normalized spatial coordinate to radius of magnetic coil
nu_ab0  = @(a,b) b.n * e^4 * a.Z^2 * b.Z^2 * ln_delt / (2*pi*eps0^2*a.m^2*b.vth^3); %scaling coefficient
nu_s    = @(v,a,b) nu_ab0(a,b) .* (1+a.m/b.m) .* psi(x(v,b.vth)) ./ x(v,b.vth); %slowing down frequency
nu_par  = @(v,a,b) nu_ab0(a,b).*(psi(x(v,b.vth))./(x(v,b.vth).^3)); %parallel diffusion frequency
nu_D    = @(v,a,b) nu_ab0(a,b).*(phi_f(x(v,b.vth)) - psi(x(v,b.vth)))./(x(v,b.vth).^3); %deflection frequency in s^-1
maxwell = @(v,offset,vth) a.n/(pi^3/2.*vth^3).*exp(-((v-offset)/vth).^2);
gauss   = @(v,x) a.n/(sqrt(2*pi)*x)*exp(-0.5*((v - x)/x).^2);
B_func = @(s) exp(s); % @(xi) B_o./(1 + xi.^2).^(3/2) %magnetic field as a function of spatial coordinate
dB_ds = @(s) exp(s); % @(xi, R_mag) -3*B_o.*xi./(R_mag.*(1 + xi.^2).^(5/2)) derivative of magnetic field
advec_time_1D = @(t) exp(-2*v_test*cos(z_test)*t);
uniform = @(x,p,t) x.*0 + 1; %uniform condition if needed

init_cond_v = @(v) maxwell(v, a.vth, 10^6);
init_cond_z = @(z) a.n*cos(z);
init_cond_s = @(s) exp(s);
init_cond_t = @(t) t*0 + 1;

analytic_solution_s = @soln_s;
    function ret = soln_s(s,p,t)
	ret = exp(s);
        if isfield(p,'norm_fac')
             ret = ret .* p.norm_fac;
        end
    end
analytic_solution_v = @soln_v;
    function ret = soln_v(v,p,t)
        ret = a.n/(pi^3/2.*v_th(b.T_eV,a.m).^3).*exp(-(v./v_th(b.T_eV,a.m)).^2);
        if isfield(p,'norm_fac')
            ret = ret .* p.norm_fac;
        end
    end

analytic_solution_z = @soln_z;
    function ret = soln_z(z,p,t)
        ret = a.n.*cos(z).*exp(-nu_D(b.vth,a,b).*t);
        if isfield(p,'norm_fac')
            ret = ret .* p.norm_fac;
        end
    end

    function ret = phi(x)
        ret = erf(x);
    end

    function ret = phi_f(x)
        ret = (x + 1./(2*x)).*erf(x) + exp(-x.^2)./sqrt(pi);
    end

    function ret = psi(x)
        dphi_dx = 2./sqrt(pi) * exp(-x.^2);
        ret = 1./(2*x.^2) .* (phi(x) - x.*dphi_dx);
        ix = find(abs(x)<1e-5); % catch singularity at boundary
        ret(ix) = 0;
    end

% throw all these into a parameters structure for use elsewhere

save_filename = 'mirror_parameters.mat';
save(save_filename);
p = load(save_filename);
delete(save_filename);

end
