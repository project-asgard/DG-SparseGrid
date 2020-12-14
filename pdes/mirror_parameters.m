function p = mirror_parameters()

% Common constants and functions for the mirror set of PDEs

% Define coil geometry:
% =========================================================================
Lx = 6;                            % [m]
Lx_offset = 0;                     % [m]
N = 200;                           % number of elements on profile
% =========================================================================
% Axial domain:
z = linspace(-Lx/2,Lx/2,N)' + Lx_offset;

m_e = 9.109*10^-31; %electron mass in kg
m_D = 3.3443*10^-27; %deuterium mass in kgs
m_H = 1.6726*10^-27; %hydrogen mass in kgs
k_b = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
e = 1.602*10^-19; %charge in Coulombs
ln_delt = 10; %Coulomb logarithm
v_th = @(T_eV,m) sqrt(2*T_eV * e/m);
eps0 = 8.85*10^-12; %permittivity of free space in Farad/m
mu0 = 4*pi*10^-7; %magnetic permeability in Henries/m
n_turns = 10;
R_mag = 0.6; %radius of individual loop in meters
I_mag = 10; %current going through individual loop in meters
B_o = mu0*n_turns*I_mag/(2*R_mag); %magnetic field under current loop
v_test = 500; %test value for velocity in coefficient
z_test = pi/2 - 1e-6; %test value for pitch angle in coefficient
advec_space_1D = v_test*cos(z_test); %advection coefficient for 1D in space
n  = [+5, +5]*1e3; %number of loops in each coil
coil_radius  = [R_mag, R_mag]; %radii of coils
z0 = [-1.5, 1.5]; %location of coils in physical space
I  = [+1.5,1.5]*1e2; %current going through each coil in amps

% species b: electrons in background
b.n = 4e19;
b.E_eV = 0;
b.T_eV = 4;
b.Z = -1;
b.m = m_e;
b.vth = v_th(b.T_eV,b.m);

%species b2: deuterium in background
b2.n = 4e19;
b2.E_eV = 0;
b2.T_eV = 4;
b2.Z = 1;
b2.m = m_D;
b2.vth = v_th(b2.T_eV,b2.m);

% species a: electrons in beam
a.n = 4e19;
a.E_eV = 1e3;
a.T_eV = 50;
a.Z = -1;
a.z0 = pi/4; %location of beam injection in pitch 
a.dz0 = a.z0/4;
a.s0 = 0;
a.ds0 = 0.1;
a.m = m_e;
a.vth = v_th(a.T_eV,a.m);
a.v_beam = sqrt(2*e.*a.E_eV./a.m);

x       = @(v,vth) v./vth; %normalized velocity to thermal velocity
xi      = @(s,R_mag) s/R_mag; %normalized spatial coordinate to radius of magnetic coil
nu_ab0  = @(a,b) b.n * e^4 * a.Z^2 * b.Z^2 * ln_delt / (2*pi*eps0^2*a.m^2*b.vth^3); %scaling coefficient
nu_s    = @(v,a,b) nu_ab0(a,b) .* (1+a.m/b.m) .* psi(x(v,b.vth)) ./ x(v,b.vth); %slowing down frequency
nu_par  = @(v,a,b) nu_ab0(a,b).*(psi(x(v,b.vth))./(x(v,b.vth).^3)); %parallel diffusion frequency
nu_D    = @(v,a,b) nu_ab0(a,b).*(phi_f(x(v,b.vth)) - psi(x(v,b.vth)))./(x(v,b.vth).^3); %deflection frequency in s^-1
maxwell = @(v,offset,vth) a.n/(pi^3/2.*vth^3).*exp(-((v-offset)/vth).^2);
gauss   = @(v,x,y) 1/(sqrt(2*pi)*x)*exp(-0.5*((v - x)/y).^2);
B_func = @(s) exp(s); % %magnetic field as a function of spatial coordinate
dB_ds = @(s) exp(s); 
% Vacuum magnetic field function:
% Function based on simple current loops:
B_func2 = @(s) sum((mu0.*n.*I./(2*coil_radius)).*(1 + ((s-z0)./coil_radius).^2 ).^(-3/2));%magnetic field from sum of loops around mirror
% Calculate normalizied vacuum magnetic field profile:
% Set the profile equal to 1 at the center of mirror:
%for ii = 1:numel(z)
%    Bz(ii) = B_func2(z(ii));
%end
dB_ds2 = @(s) sum(-3*(mu0.*n.*I./(2*coil_radius)).*(s-z0)./coil_radius./(coil_radius.*(1 + ((s-z0)./coil_radius).^2).^(5/2))); %derivative of field around loops
for ii = 1:numel(z)
    dBdz(ii) = dB_ds2(z(ii));
end
advec_time_1D = @(t) exp(-2*v_test*cos(z_test)*t);
uniform = @(x,p,t) x.*0 + 1; %uniform condition if needed

init_cond_v = @(v,p,t) maxwell(v,a.v_beam,a.vth);
init_cond_z = @(z,p,t) gauss(z,a.z0,a.dz0);
init_cond_s = @(s,p,t) gauss(s,a.s0,a.ds0);
init_cond_t = @(t,p) t*0 + 1;

boundary_cond_v = @(v,p,t) maxwell(v,0,b.vth);%exp(-nu_D(v,a,b).*t);
boundary_cond_z = @(z,p,t) z.*0 + 0.153667;
boundary_cond_s = @(s,p,t) s.*0 + 1;
boundary_cond_t = @(t,p) t.*0 + 1;

s_z = @(z,p,t) -0.5.*(sin(z)+cot(z).*cos(z)); %source pitch fuction for mirror3
s_v = @(v,p,t) nu_D(v,a,b);
s_space = @(s,p,t) s.*0 + 1;
s_time = @(t,p) t.*0 + 1;

source_3D = {s_z, s_v, s_space, s_time};


soln_s = @solution_s;
    function ret = solution_s(s,p,t)
	ret = exp(s);
    end
soln_v = @solution_v;
    function ret = solution_v(v,p,t)
        ret = a.n/(pi^3/2.*v_th(b.E_eV,a.m).^3).*exp(-(v./v_th(b.E_eV,a.m)).^2);
    end
soln_z = @solution_z;
    function ret = solution_z(z,p,t)
        ret = cos(z/2);
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
