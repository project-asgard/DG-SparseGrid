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

x       = @(v,vth) v./vth;
nu_ab0  = @(a,b) b.n * e^4 * a.Z^2 * b.Z^2 * ln_delt / (2*pi*eps0^2*a.m^2*b.vth^3); %scaling coefficient
nu_s    = @(v,a,b) nu_ab0(a,b) .* (1+a.m/b.m) .* psi(x(v,b.vth)) ./ x(v,b.vth); %slowing down frequency
nu_par  = @(v,a,b) nu_ab0(a,b).*(psi(x(v,b.vth))./(x(v,b.vth).^3)); %parallel diffusion frequency
nu_D    = @(v,a,b) nu_ab0(a,b).*(phi_f(x(v,b.vth)) - psi(x(v,b.vth)))./(x(v,b.vth).^3); %deflection frequency in s^-1
maxwell = @(v,offset,vth) a.n/(pi^3/2.*vth^3).*exp(-((v-offset)/vth).^2);
gauss   = @(v,x) a.n/(sqrt(2*pi)*x)*exp(-0.5*((v - x)/x).^2);

init_cond_v = @(v) maxwell(v,1.5e7,v_th(a.T_eV,a.m));
init_cond_z = @(z) z.*0 + 1;
init_cond_t = @(t) t*0 + 1;

analytic_solution_v = @soln_v;
function ret = soln_v(v,p,t)
ret = a.n/(pi^3/2.*v_th(b.T_eV,a.m).^3).*exp(-(v./v_th(b.T_eV,a.m)).^2);
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

end