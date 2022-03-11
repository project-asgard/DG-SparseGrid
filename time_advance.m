function [f,sol] = time_advance(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax)

sol = 0;

if strcmp(opts.timestep_method,'BE')
    % Backward Euler (BE) first order
    f = backward_euler(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'FE')
    % Forward Euler (FE) first order
    f = forward_euler(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'matrix_exponential')
    % Matrix Exponential (ME) all order
    f = matrix_exponential(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'time_independent')
    % time independent d/dt==0
    f = time_independent(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'CN')
    % Crank Nicolson (CN) second order
    f = crank_nicolson(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'IMEX')
    % IMEX
    f = imex(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif sum(strcmp(opts.timestep_method,{'ode15i','ode15s','ode45','ode23s'}))>0
    % One of the atlab ODE integrators.
    [f,sol] = ODEm(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
else
    % RK3 TVD
    f = RungeKutta3(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
end

end

%% Use Matlab ODE solvers

function [fout,sol] = ODEm(pde,opts,A_data,f0,t0,dt,deg,hash_table,Vmax,Emax)

assert(opts.time_independent_A,...
    'Since get_coeff_mats is not called within the ode RHS function this presently only works for time_independent_A,true');

% when using the matlab ODE integrators we just pass them the start and
% final time (when 'time_independent_A','true') and have it return the
% intermeadiate times

tf = (opts.num_steps+1) * dt; % the +1 here is to account for very small deviations in the end time requested from sol
times = [t0:dt:tf];

clear strCR;

stats = 'off';
output_func = '';
if(~opts.quiet)
    stats = 'on';
    output_func = @odetpbar;
end

applyLHS = ~isempty(pde.termsLHS);

if applyLHS
    [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg,Vmax,Emax);
else
    [~,A] = apply_A(pde,opts,A_data,f0,deg,Vmax,Emax);
end

options = odeset('RelTol',1e-6,'AbsTol',1e-8,'Stats',stats,'OutputFcn',output_func,'Jacobian',A);

    function dfdt = explicit_ode(t,f)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
        if applyLHS
%             [f1,~,ALHS] = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
            f1 = ALHS\(A*f + source + bc);
            dfdt = ALHS\(f1 + source + bc);       
        else
%             [f1] = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
            f1 = A*f;
            dfdt = f1 + source + bc;       
        end
    end

    function res = implicit_ode(t,f,dfdt)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
%         f1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
        f1 = A*f;
        res = dfdt - (f1 + source + bc);
    end

if strcmp(opts.timestep_method,'ode45')
    
    sol = ode45(@explicit_ode,times,f0,options);
    fout = deval(sol,t0+dt);
    
elseif strcmp(opts.timestep_method,'ode15s')
    
     sol = ode15s(@explicit_ode,times,f0,options);
     fout = deval(sol,t0+dt);
    
elseif strcmp(opts.timestep_method,'ode23s')
    
    sol = ode23s(@explicit_ode,times,f0,options);
    fout = deval(sol,t0+dt);
    
elseif strcmp(opts.timestep_method,'ode15i')

    if applyLHS
        error('ERROR: ode15i not yet implemented for LHS=true');
    end
    
    dfdt0 = f0.*0;
    [f0,dfdt0,resnrm] = decic(@implicit_ode,t0,f0,f0.*0+1,dfdt0,[]);
    sol = ode15i(@implicit_ode,times,f0,dfdt0);
    fout = deval(sol,t0+dt);
end

end

%% 3-rd Order Kutta Method (explicit time advance)
function fval = RungeKutta3(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax)

dims = pde.dimensions;

%%
% Sources
c2 = 1/2; c3 = 1;
source1 = source_vector(pde,opts,hash_table,t);
source2 = source_vector(pde,opts,hash_table,t+c2*dt);
source3 = source_vector(pde,opts,hash_table,t+c3*dt);

%%
% Inhomogeneous dirichlet boundary conditions
bc1 = boundary_condition_vector(pde,opts,hash_table,t);
bc2 = boundary_condition_vector(pde,opts,hash_table,t+c2*dt);
bc3 = boundary_condition_vector(pde,opts,hash_table,t+c3*dt);

% %%
% Apply any non-identity LHS mass matrix coefficient

applyLHS = ~isempty(pde.termsLHS);

a21 = 1/2; a31 = -1; a32 = 2;
b1 = 1/6; b2 = 2/3; b3 = 1/6;

if applyLHS
    [k1,A1,ALHS] = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
    rhs1 = source1 + bc1;
    %     invMatLHS = inv(ALHS); % NOTE : assume time independent for now for speed.
    %     k1 = invMatLHS * (k1 + rhs1);
    k1 = ALHS \ (k1 + rhs1);
    y1 = f + dt*a21*k1;
    
    [k2] = apply_A(pde,opts,A_data,y1,deg,Vmax,Emax);
    rhs2 = source2 + bc2;
    %     k2 = invMatLHS * (k2 + rhs2);
    k2 = ALHS \ (k2 + rhs2);
    y2 = f + dt*(a31*k1 + a32*k2);
    
    k3 = apply_A(pde,opts,A_data,y2,deg,Vmax,Emax);
    rhs3 = source3 + bc3;
    %     k3 = invMatLHS * (k3 + rhs3);
    k3 = ALHS \ (k3 + rhs3);
else
    k1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax)  + source1 + bc1;
    y1 = f + dt*a21*k1;
    k2 = apply_A(pde,opts,A_data,y1,deg,Vmax,Emax) + source2 + bc2;
    y2 = f + dt*(a31*k1 + a32*k2);
    k3 = apply_A(pde,opts,A_data,y2,deg,Vmax,Emax) + source3 + bc3;
end

fval = f + dt*(b1*k1 + b2*k2 + b3*k3);

end

%% Time independent solve d/dt==0
function f1 = time_independent(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

[~,A] = apply_A(pde,opts,A_data,f0,deg);

f1 = -A \ (s0+bc0);

end

%% Matrix Exponential
function f1 = matrix_exponential(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

if applyLHS
    [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
    f1 = expm(dt*(ALHS\A))*f0;
else
    [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
    f1 = expm(A*dt)*f0;
end

Q = A*f0;


function ret = rhs_integrand(u)
    s = source_vector(pde,opts,hash_table,u);
    bc = boundary_condition_vector(pde,opts,hash_table,u);
    uu = t+dt-u;
    if applyLHS
        ret = ALHS\(expm(uu*(ALHS\A)) * (s+bc));
    else
        ret = expm(uu*A) * (s+bc);
    end
end

rhs_integral = integral(@rhs_integrand,t,t+dt,'ArrayValued',true);%,'RelTol',1e-5,'AbsTol',1e-5);
f1 = f1 + rhs_integral;

end

%% Forward Euler
function f1 = forward_euler(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

applyLHS = ~isempty(pde.termsLHS);

if applyLHS
    error('apply LHS not implemented for FE');
else
    f1 = f0 + dt * (apply_A(pde,opts,A_data,f0,deg) + s0 + bc0);
end

end

%% Backward Euler (first order implicit time advance)
function f1 = backward_euler(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

build_A = true;
persistent dA dALHS A ALHS;

if opts.time_independent_A
    if ~isempty(dA)
        build_A = false;
    end
end
if opts.time_independent_build_A
    if ~isempty(A)
        build_A = false;
    end
end

if build_A
    if applyLHS
        [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);       
    else
        [~,A] = apply_A(pde,opts,A_data,f0,deg);
    end
else
end

if ~opts.time_independent_A || isempty(dA)
    if applyLHS
        AA = ALHS - dt*A;
    else
        I = speye(numel(diag(A)));
        AA = I - dt*A;
    end
end

if applyLHS
    b  = ALHS * (f0 + dt*(s0 + bc0));  
else
    b  = f0 + dt*(s0 + bc0);   
end

if opts.time_independent_A
    if isempty(dA)
        if ~opts.quiet;disp(['   Precomputing matrix decomposition ...']);end
        if ~opts.quiet;disp(['     nnz/numel*100 (before droptol) = ', num2str(nnz(AA)/numel(AA)*100,'%2.0f'),'%']);end
        AA(abs(AA)<1e-16)=0;
        if ~opts.quiet;disp(['     nnz/numel*100  (after droptol) = ', num2str(nnz(AA)/numel(AA)*100,'%2.0f'),'%']);end
        dA = decomposition(AA);
        clear A AA;
        if ~opts.quiet;disp('    DONE');end
    end
end

if opts.time_independent_A
    f1 = dA \ b;
else
    f1 = AA \ b;
end
 
end


%% Crank Nicolson (second order implicit time advance)
function f1 = crank_nicolson(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

s0 = source_vector(pde,opts,hash_table,t);
s1 = source_vector(pde,opts,hash_table,t+dt);

bc0 = boundary_condition_vector(pde,opts,hash_table,t);
bc1 = boundary_condition_vector(pde,opts,hash_table,t+dt);

build_A = true;
persistent dA dALHS A ALHS;

if opts.time_independent_A
    if ~isempty(dA)
        build_A = false;
    end
end
if opts.time_independent_build_A
    if ~isempty(A)
        build_A = false;
    end
end

if build_A
    if applyLHS
        [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
    else
        [~,A] = apply_A(pde,opts,A_data,f0,deg);      
    end
end

if ~opts.time_independent_A || isempty(dA)
if applyLHS
    AA = 2*ALHS - dt*A;
else
    I = speye(numel(diag(A)));
    AA = 2*I - dt*A;
end
end

if applyLHS
    b  = ALHS * (2*f0 + dt*A*f0 + dt*(s0+s1+bc0+bc1));  
else
    b  = 2*f0 + dt*A*f0 + dt*(s0+s1+bc0+bc1);  
end

if opts.time_independent_A
    if isempty(dA)
        if ~opts.quiet;disp(['   Precomputing matrix decomposition ...']);end
        if ~opts.quiet;disp(['     nnz/numel*100 (before droptol) = ', num2str(nnz(AA)/numel(AA)*100,'%2.0f'),'%']);end
        AA(abs(AA)<1e-16)=0;
        if ~opts.quiet;disp(['     nnz/numel*100  (after droptol) = ', num2str(nnz(AA)/numel(AA)*100,'%2.0f'),'%']);end
        dA = decomposition(AA);
        clear AA;
        if ~opts.quiet;disp('    DONE');end
    end
end

if opts.time_independent_A
    f1 = dA \ b;
else
    f1 = AA \ b;
end

end

function f1 = imex(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

persistent pde_1d
persistent nodes
persistent Meval

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% C++ Implementation Notes %%%%%%%%%%%%%%%%%%%%%%%

% the fast_2d_matrix_apply routine should not be 
% implemented in the C++ version.  It is only there to
% make the MATLAB version viable for 1+1 runs.

% Use the standard apply_A routine instead.  It will need
% to be modified to incorporate the imex_flags.  This is
% already done in the matlab version.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%Switch for IMEX iteration: 
% 1 : Backward Euler in explicit terms 'E',
%     Backward Euler in implicit terms 'I'.
% 0 : SSP-RK2 IMEX.  'E' terms are actually handled explicitly.

%%% C++ : Focus on the case BEFE = 0

BEFE = 0;

if isempty(Meval)
    %Compute everything that will not change per time iteration

    %Get quadrature points in realspace stiffness matrix calculation
    [Meval,nodes] = matrix_plot_D(pde,opts,pde.dimensions{1});
    
    %Make fake pde file for use in physical realspace transformation
    pde_1d.dimensions = pde.dimensions(1);
        
end

%Create moment matrices that take DG function in (x,v) and transfer it
%to DG function in x.
if numel(pde.dimensions) >= 2
    moment_mat = cell(numel(pde.moments),1);
    for i=1:numel(pde.moments)
        moment_mat{i} = moment_reduced_matrix(opts,pde,A_data,hash_table,i);
    end
end

hash_table_1D = hash_table_2D_to_1D(hash_table,opts);

assert(isempty(pde.termsLHS),'LHS terms currently not supported by IMEX');
%Ignoring sources for now
assert(isempty(pde.sources),'Sources currently not supported by IMEX');

if BEFE

    %Explicit Update
    [f_E,flag,relres,iter] = bicgstabl(@(x) x - dt*fast_2d_matrix_apply(opts,pde,A_data,x,'E'),f0,1e-12,numel(f0));
    if flag ~= 0
        fprintf('BICGSTABL did not converge.  flag = %d, relres = %5.4e\n',flag,relres);
        assert(relres < 1e-10)
    end

    %Now get moments
    mom0 = moment_mat{1}*f_E; %integral of (f,1)_v
    mom0_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom0,hash_table_1D);
    pde.params.n  = @(x) interp1(nodes,mom0_real,x,'linear','extrap');

    mom1 = moment_mat{2}*f_E; %integral of (f,v)_v
    mom1_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom1,hash_table_1D);
    pde.params.u  = @(x) interp1(nodes,mom1_real,x,'linear','extrap')./pde.params.n(x);

    mom2 = moment_mat{3}*f_E; %integral of (f,v^2)_v
    mom2_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom2,hash_table_1D);
    pde.params.th = @(x) interp1(nodes,mom2_real,x,'linear','extrap')./pde.params.n(x) - pde.params.u(x).^2;

    if ~opts.quiet
        figure(1000);
        subplot(2,2,1);
        plot(nodes,pde.params.n(nodes));
        title('n_f');
        subplot(2,2,2);
        plot(nodes,pde.params.u(nodes));
        title('u_f');
        subplot(2,2,3);
        plot(nodes,pde.params.th(nodes));
        title('th_f');
        sgtitle("Fluid Variables. t = "+num2str(t+dt));
    end

    %Update coefficients
    pde = get_coeff_mats(pde,opts,t,0);

    %Implicit Update
    [f1,flag,relres,iter] = bicgstabl(@(x) x - dt*fast_2d_matrix_apply(opts,pde,A_data,x,'I'),f_E,1e-12,numel(f0));
    if flag ~= 0
        fprintf('BICGSTABL did not converge.  flag = %d, relres = %5.4e\n',flag,relres);
        assert(relres < 1e-10)
    end

else %%Trying imex deg 2 version
    %Here f1 = f^{n+1}, f0 = f^n
    
    %%%%%
    %%% First stage
    %%%%%
    
    % nothing happens  f_1 = f_n
    
    %%%%%
    %%% Second stage
    %%%%%
    
    %Explicit step
    f_2s = f0 + dt*fast_2d_matrix_apply(opts,pde,A_data,f0,'E');
    
    %Create rho_2s
    mom0 = moment_mat{1}*f_2s; %integral of (f,1)_v
    mom0_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom0,hash_table_1D);
    pde.params.n  = @(x) interp1(nodes,mom0_real,x,'nearest','extrap');

    mom1 = moment_mat{2}*f_2s; %integral of (f,v)_v
    mom1_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom1,hash_table_1D);
    pde.params.u  = @(x) interp1(nodes,mom1_real,x,'nearest','extrap')./pde.params.n(x);

    mom2 = moment_mat{3}*f_2s; %integral of (f,v^2)_v
    mom2_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom2,hash_table_1D);
    pde.params.th = @(x) interp1(nodes,mom2_real,x,'nearest','extrap')./pde.params.n(x) - pde.params.u(x).^2;
        
    %Update coefficients
    pde = get_coeff_mats(pde,opts,t,0);
    
    %f2 now
    [f_2,flag,relres,iter] = bicgstabl(@(x) x - dt*fast_2d_matrix_apply(opts,pde,A_data,x,'I'),f_2s,1e-12,numel(f0));
    if flag ~= 0
        fprintf('BICGSTABL did not converge.  flag = %d, relres = %5.4e\n',flag,relres);
        assert(relres < 1e-10)
    end
    
    %%%%%
    %%% Third stage
    %%%%%
    
    f_3s = f0 + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f0+f_2,'E') ...
              + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f_2,'I');
    
    %Create rho_3s
    mom0 = moment_mat{1}*f_3s; %integral of (f,1)_v
    mom0_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom0,hash_table_1D);
    %mom0_vals = singleD_to_multiD(num_dims-1,mom0_real,nodes(1));
    pde.params.n  = @(x) interp1(nodes,mom0_real,x,'nearest','extrap');

    mom1 = moment_mat{2}*f_3s; %integral of (f,v)_v
    mom1_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom1,hash_table_1D);
    pde.params.u  = @(x) interp1(nodes,mom1_real,x,'nearest','extrap')./pde.params.n(x);

    mom2 = moment_mat{3}*f_3s; %integral of (f,v^2)_v
    mom2_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom2,hash_table_1D);
    pde.params.th = @(x) interp1(nodes,mom2_real,x,'nearest','extrap')./pde.params.n(x) - pde.params.u(x).^2;
        
    %Update coefficients
    pde = get_coeff_mats(pde,opts,t,0);

    [f_3,flag,relres,iter] = bicgstabl(@(x) x - dt*fast_2d_matrix_apply(opts,pde,A_data,x,'I'),f_3s,1e-12,numel(f0));
    if flag ~= 0
        fprintf('BICGSTABL did not converge.  flag = %d, relres = %5.4e\n',flag,relres);
        assert(relres < 1e-10)
    end
    
    %Update timestep to final stage
    f1 = f_3;
    
    %Plot moments
    if ~opts.quiet || 1
        
        fig1 = figure(1000);
        fig1.Units = 'Normalized';
        fig1.Position = [0.5 0.5 0.3 0.3];
        subplot(2,2,1);
        plot(nodes,pde.params.n(nodes));
        title('n_f');
        subplot(2,2,2);
        plot(nodes,pde.params.u(nodes));
        title('u_f');
        subplot(2,2,3);
        plot(nodes,pde.params.th(nodes));
        title('th_f');
        sgtitle("Fluid Variables. t = "+num2str(t+dt));
        drawnow
        
    end
    
end

end
