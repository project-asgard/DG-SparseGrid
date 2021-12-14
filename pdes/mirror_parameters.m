function p = mirror_parameters(opts)

% Common constants and functions for the mirror set of PDEs

% Define coil geometry:
% =========================================================================
Lx = 6;                            % [m]
Lx_offset = 0;                     % [m]
N = 200;                           % number of elements on profile
% =========================================================================
% Axial domain:
%z = linspace(-Lx/2,Lx/2,N)' + Lx_offset;
%speed of light

c = 3*10^8; %m/s

%setting up integrand matrix

gamma = @(u) sqrt(1 + u.^2./c^2); %relativistic correction

m_e = 9.109*10^-31; %electron mass in kg
m_D = 3.3443*10^-27; %deuterium mass in kgs
m_H = 1.6726*10^-27; %hydrogen mass in kgs
boltz = 1.380*10^-23; %Boltzmann constant in Joules/Kelvin
e = 1.602*10^-19; %charge in Coulombs
ln_delt = 10; %Coulomb logarithm
v_th = @(T_eV,m) sqrt(2*T_eV * e/m);
eps0 = 8.85*10^-12; %permittivity of free space in Farad/m
mu0 = 4*pi*10^-7; %magnetic permeability in Henries/m
n_turns = 10;
radius_loop = 0.6; %radius of individual loop in meters
current_loop = 10; %current going through individual loop in meters
magn_field = mu0*n_turns*current_loop/(2*radius_loop); %magnetic field under current loop
vel_test = 500; %test value for velocity in coefficient
pitch_test = pi/2 - 1e-6; %test value for pitch angle in coefficient
advec_space_1D = vel_test*cos(pitch_test); %advection coefficient for 1D in space
alpha_z = @(z) z.*0;
coil_nloops  = [+5, +5]*1e3; %number of loops in each coil
coil_radius  = [radius_loop, radius_loop]; %radii of coils
coil_coords = [-1.5, 1.5]; %location of coils in physical space
coil_i  = [+1.5,1.5]*1e2; %current going through each coil in amps


vel_norm       = @(v,vth) v./vth; %normalized velocity to thermal velocity
space_norm      = @(s,R_mag) s/R_mag; %normalized spatial coordinate to radius of magnetic coil
nu_ab0  = @(a,b) b.n * e^4 * a.Z^2 * b.Z^2 * ln_delt / (2*pi*eps0^2*a.m^2*b.vth^3); %scaling coefficient
nu_s    = @(v,a,b) nu_ab0(a,b) .* (1+a.m/b.m) .* psi(vel_norm(v,b.vth)) ./ vel_norm(v,b.vth); %slowing down frequency
nu_par  = @(v,a,b) nu_ab0(a,b).*(psi(vel_norm(v,b.vth))./(vel_norm(v,b.vth).^3)); %parallel diffusion frequency
nu_D    = @(v,a,b) nu_ab0(a,b).*(phi(vel_norm(v,b.vth)) - psi(vel_norm(v,b.vth)))./(vel_norm(v,b.vth).^3); %deflection frequency in s^-1
%nu_E    = @(v,a,b) nu_ab0(a,b).*(2.*(a.m/b.m).*psi(vel_norm(v,b.vth))./(vel_norm(v,b.vth)) - dphidx(vel_norm(v,b.vth))./(vel_norm(v,b.vth)).^2);
nu_E    = @(v,a,b) nu_ab0(a,b).*(1./vel_norm(v,b.vth).^3).*(a.m/b.m).*(phi(vel_norm(v,b.vth)) - (1 + b.m/a.m).*vel_norm(v,b.vth).*dphidx(vel_norm(v,b.vth))); %Hinton definition of nu_E
maxwell = @(v,offset,vth,a) a.n/(pi^(3/2).*vth^3).*exp(-((v-offset)/vth).^2);
gauss   = @(v,x,y,a) a.n/(sqrt(2*pi)*y)*exp(-0.5*((v - x)/y).^2);
B_func = @(s) exp(s); % %magnetic field as a function of spatial coordinate
dB_ds = @(s) exp(s); 
E_dreicer_si = 1;

switch opts.case_
	case 1
        % species b: electrons in background
        b.n = 4e19; %number density of background species in m^-3
        b.E_eV = 0;%energy of background species
        b.T_eV = 50;%temperature/spread of background species
        b.Z = -1;%atomic number of species
        b.m = m_e; %mass of background
        b.vth = v_th(b.T_eV,b.m);%thermal velocity of background
        
        %species b2: deuterium in background
        b2.n = 4e19;
        b2.E_eV = 0;
        b2.T_eV = 50;
        b2.Z = 1;
        b2.m = m_D;
        b2.vth = v_th(b2.T_eV,b2.m);
        
        % species a: species in beam
        a.n = 4e19;
        a.E_eV = 3e3; %energy of beam species
        a.T_eV = 50; %temperature/spread of beam species
        a.Z = 1;
        a.z0 = pi/4; %location of beam injection in pitch angle (radians)
        a.dz0 = sqrt(a.T_eV/a.E_eV); %spread of beam in pitch angle
        a.s0 = 0; %location of beam in space (m)
        a.ds0 = 2.5; %spread of beam in space (m)
        a.m = m_D;
        a.vth = v_th(a.T_eV,a.m);
a.v_beam = sqrt(2*e.*a.E_eV./a.m); %velocity of beam species
	case 2
 		E = 10^-4*E_dreicer_si;
	case 3
        n_cgs = 8e14; %equilibrium density in cm.^-3
        m_e_cgs = 9.109*10^-28; %electron mass in g
        m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
        m_He_cgs = 6.7*10^-24; %helium 4 mass in g
        m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
        m_Ne_cgs = 3.3509177*10^-23; %Neon mass in g
        temp_cgs = 1.6022e-10; %temperature in erg
        ln_delt = 18; %Coulomb logarithm
        %         params_cgs.a.vth = sqrt(2*temp_cgs/m_e_cgs);
        %         params_cgs.b.vth = sqrt(2*temp_cgs/m_e_cgs);
        params_cgs.a.m = m_e_cgs; %beam is electrons
        params_cgs.b.m = m_Ne_cgs; %background ions
        params_cgs.b2.m = m_e_cgs; %background electrons
        params_cgs.a.Z = -1;
        params_cgs.b.Z = 1;
        params_cgs.b2.Z = -1;
        params_cgs.a.n = n_cgs;
        params_cgs.b.n = n_cgs/params_cgs.b.Z;
        params_cgs.b2.n = n_cgs;
        params_cgs.e = 4.803*10^-10; %charge in Fr
        params_cgs.E = 2.6e-7; %E field in statvolt/cm
        params_cgs.v_c = @(E_cgs) sqrt(4*pi*params_cgs.e^3 * params_cgs.a.n * ln_delt / (params_cgs.a.m * E_cgs));
        a.n = 10^6*params_cgs.a.n;%converting to m^-3
        b.n = 10^6*params_cgs.b.n;
        b2.n = 10^6*params_cgs.b2.n;
        %         params.a.vth = 0.01*params_cgs.a.vth; %converting to m/s
        %         params.b.vth = 0.01*params_cgs.b.vth;
        a.m = 0.001*params_cgs.a.m; %converting to kg
        b.m = 0.001*params_cgs.b.m;
        b2.m = 0.001*params_cgs.b2.m;
        a.Z = params_cgs.a.Z;
        b.Z = params_cgs.b.Z;
        b2.Z = params_cgs.b2.Z;
        %E = 2.9979*10^4*params_cgs.E; %converting to V/m
        a.E_eV = 1000;
        a.T_eV = 5.11*10^3;
        b.T_eV = a.T_eV;
        b2.T_eV = a.T_eV;
        a.vth = v_th(a.T_eV,a.m);
        a.v_beam = 0;%a.vth/6;
        b2.vth = v_th(b2.T_eV,b2.m);
        E_dreicer_si = a.n.*e^3*ln_delt/(2*pi*eps0^2*a.m ...
            *a.vth^2);
        frac = 0.05;
        E = frac*E_dreicer_si;
        %vel_norm = @(v,vth) v./vth; %normalized velocity to thermal velocity
        maxwell = @(v,offset,vth,a) a.n/(pi.^(3/2)*vth^3).*exp(-((v-offset)/vth).^2);
        soln_v = @(v,p,t) solution_v(v,p,t);
        f0_v = @(v) maxwell(v,a.v_beam,a.vth,a);
        init_cond_v = @(v,p,t) f0_v(v);
        %params_cgs.nu_ab0  = @(a,b) b.n * params_cgs.e^4 * a.Z^2 * b.Z^2 * params_cgs.ln_delt / (pi^3/2.*a.m^2*b.vth^3); %scaling coefficient
        %params.eps0 = 1/(4*pi);
        if isfield(opts.cmd_args,'delta')
            delta = opts.cmd_args.delta;
        end
        if isfield(opts.cmd_args,'E')
            E = opts.cmd_args.E;
        end
        if isfield(opts.cmd_args,'frac')
            frac = opts.cmd_args.frac;
        end
        
        if isfield(opts.cmd_args,'Z')
            b.Z = opts.cmd_args.Z;
        end
        if isfield(opts.cmd_args,'tau')
           tau = opts.cmd_args.tau;
        end
        if isfield(opts.cmd_args,'nu_ee')
            nu_ee = opts.cmd_args.nu_ee;
        end      
        if isfield(opts.cmd_args,'m')
            b.m = opts.cmd_args.m;
        end
        if isfield(opts.cmd_args,'n')
            a.n = opts.cmd_args.n;
            b.n = opts.cmd_args.n/opts.cmd_args.Z;
            b2.n = opts.cmd_args.n;
        end
        E = frac*E_dreicer_si;
        b.vth = v_th(b.T_eV,b.m);
        disp(['E: ', num2str(E)]);
        disp(['Z: ', num2str(b.Z)]);
        disp(['E/E_D:', num2str(frac)]);
    case 4
        n_cgs = 5e13; %equilibrium density in cm.^-3
        m_e_cgs = 9.109*10^-28; %electron mass in g
        m_D_cgs = 3.3443*10^-24; %Deuterium mass in g
        m_He_cgs = 6.7*10^-24; %helium 4 mass in g
        m_B_cgs = 1.82*10^-23; %Boron 11 mass in g
        temp_cgs = 1.6022e-10; %temperature in erg
        params_cgs.a.n = n_cgs;
        params_cgs.b.n = n_cgs;
        params_cgs.b2.n = n_cgs;
        params_cgs.a.m = m_e_cgs; %beam is electrons
        params_cgs.b.m = m_e_cgs; %background is electrons
        params_cgs.b2.m = m_B_cgs;
        params_cgs.a.Z = -1;
        params_cgs.b.Z = -1;
        params_cgs.b2.Z = 1;
        params_cgs.e = 4.803*10^-10; %charge in Fr
        params_cgs.E = 2.6e-5; %E field in statvolt/cm
        a.n = 10^6*params_cgs.a.n;%converting to m^-3
        b.n = 10^6*params_cgs.b.n;
        b2.n = 10^6*params_cgs.b2.n;
        
        a.m = 0.001*params_cgs.a.m; %converting to kg
        b.m = 0.001*params_cgs.b.m;
        b2.m = 0.001*params_cgs.b2.m;
        a.Z = params_cgs.a.Z;
        b.Z = params_cgs.b.Z;
        b2.Z = params_cgs.b2.Z;
        ln_delt = 27.0857;
        a.T_eV = 1.33e4;%2/3*params.a.E_eV;
        b.T_eV = a.T_eV;
        b2.T_eV = a.T_eV;
        a.v_beam = 0;
        a.vth = v_th(a.T_eV,a.m);
        b.vth = v_th(b.T_eV,b.m);
        b2.vth = v_th(b2.T_eV,b2.m);
        E_dreicer_si = a.n.*e^3*ln_delt/(2*pi*eps0^2*a.m ...
            *a.vth^2);
        E = 10^-2*E_dreicer_si;
        %vel_norm = @(v,vth) v./vth; %normalized velocity to thermal velocity
        maxwell = @(v,offset,vth,a) a.n/(pi.^(3/2)*vth^3).*exp(-((v-offset)/vth).^2);
        init_cond_v = @(v,p,t) maxwell(v,0,a.vth,a);
        soln_v = @(v,p,t) solution_v(v,p,t);
        %params_cgs.nu_ab0  = @(a,b) b.n * params_cgs.e^4 * a.Z^2 * b.Z^2 * params_cgs.ln_delt / (pi^3/2.*a.m^2*b.vth^3); %scaling coefficient
        %params.eps0 = 1/(4*pi);
end

% Vacuum magnetic field function:
% Function based on simple current loops:

B_func2 = @get_B;

    function B = get_B(s)
        B = s.*0;
        for coil = 1:numel(coil_radius)
            zz = coil_coords(coil);
            rr = coil_radius(coil);
            II = coil_i(coil);
            nn = coil_nloops(coil);
            B = B + (mu0.*nn.*II./(2*rr)).*(1 + ((s-zz)./rr).^2 ).^(-3/2);%magnetic field from sum of loops around mirror
        end
    end

dB_ds2 = @get_dBds; 

    function dBds = get_dBds(s)
        
        dBds = s.*0;
        for coil = 1:numel(coil_radius)
            zz = coil_coords(coil);
            rr = coil_radius(coil);
            II = coil_i(coil);
            nn = coil_nloops(coil);    
            dBds = dBds -3*(mu0.*nn.*II./(2.*rr)).*(s-zz)./rr./(rr.*(1 + ((s-zz)./rr).^2).^(5/2)); %derivative of field around loops
        end
    end

advec_time_1D = @(t) exp(-2*vel_test*cos(pitch_test)*t);
uniform = @(x,p,t) x.*0 + 1; %uniiform condition if needed

f0_v = @f_init_v;
    function res = f_init_v(v)
        %res = zeros(size(v));
        res = maxwell(v,a.v_beam,a.vth,a);
        %res = gauss(v,a.v_beam,a.vth,a);
    end

init_cond_v = @(v,p,t) maxwell(v,a.v_beam,a.vth,a);
init_cond_z = @(z,p,t) z.*0 + 1;%gauss(z,a.z0,a.dz0,a);
init_cond_s = @(s,p,t) gauss(s,a.s0,a.ds0,a);
init_cond_t = @(t,p) t*0 + 1;

boundary_cond_v = @(v,p,t) maxwell(v,0,b.vth,a);%exp(-nu_D(v,a,b).*t);
boundary_cond_z = @(z,p,t) z.*0 + 1;
boundary_cond_s = @(s,p,t) s.*0;
boundary_cond_t = @(t,p) t.*0 + 1;

source_pitch = @(z,p,t) -0.5.*(sin(z)+cot(z).*cos(z)); %source pitch fuction for mirror3
source_vel = @(v,p,t) nu_D(v,a,b);
source_space = @(s,p,t) s.*0 + 1;
source_time = @(t,p) t.*0 + 1;

source_3D = {source_pitch, source_vel, source_space, source_time};


soln_s = @solution_s;
    function ret = solution_s(s,p,t)
	ret = exp(s);
    end
soln_v = @solution_v;
    function ret = solution_v(v,p,t)
%        ret = init_cond_v(v);
         ret = maxwell(v,0,a.vth,a);
%        ret = a.n/(pi^3/2.*v_th(b.T_eV,a.m).^3).*exp(-(v./v_th(b.T_eV,a.m)).^2);
        if isfield(p,'norm_fac')
            ret = p.norm_fac .* ret;
        end
    end
soln_z = @solution_z;
    function ret = solution_z(z,p,~)
        ret = z.*0 + 1;
    end

    function ret = phi(x)
        ret = erf(x);
    end

    function ret = dphidx(x)
        ret = 2./sqrt(pi) * exp(-x.^2);
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
