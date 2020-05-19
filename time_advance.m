function f = time_advance(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax)

if strcmp(opts.timestep_method,'BE')
    % Backward Euler (BE) first order
    f = backward_euler(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'CN')
    % Crank Nicolson (CN) second order
    f = crank_nicolson(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif sum(strcmp(opts.timestep_method,{'ode15i','ode15s','ode45'}))>0
    % One of the atlab ODE integrators.
    f = ODEm(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);  
else
    % RK3 TVD
    f = RungeKutta3(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
end

end

%% Use Matlab ODE solvers

function fval = ODEm(pde,opts,A_data,f0,t0,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);
if applyLHS
    error('ERROR: Matlab ODE integrators not yet implemented for LHS=true');
end

    function dfdt = explicit_ode(t,f)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
        f1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
        dfdt = f1 + source + bc;
    end

    function res = implicit_ode(t,f,dfdt)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
        f1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
        res = dfdt - (f1 + source + bc);
    end

if strcmp(opts.timestep_method,'ode45')
    
    if(~opts.quiet);disp('Using ode45');end
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'Stats','off');
    [tout,fout] = ode45(@explicit_ode,[t0 t0+dt],f0,options);
    
elseif strcmp(opts.timestep_method,'ode15s')
    
    S = eye(numel(f0));
    if(~opts.quiet);disp('Using ode15s');end
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,...
        'Stats','off','OutputFcn',@odetpbar,'Refine',20);%,'JPattern',S);
    [tout,fout] = ode15s(@explicit_ode,[t0 t0+dt],f0,options);
    
elseif strcmp(opts.timestep_method,'ode15i')
    
    dfdt0 = f0.*0;
    [f0,dfdt0,resnrm] = decic(@implicit_ode,t0,f0,f0.*0+1,dfdt0,[]);
    if(~opts.quiet);disp('Using ode15i');end
    [tout,fout] = ode15i(@implicit_ode,[t0 t0+dt],f0,dfdt0);
    
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


%% Backward Euler (first order implicit time advance)
function f1 = backward_euler(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

% %%
% Apply any non-identity LHS mass matrix coefficient

applyLHS = ~isempty(pde.termsLHS);

[~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);

I = eye(numel(diag(A)));

if applyLHS
    %     AA = I - dt*inv(ALHS)*A;
    %     b = f0 + dt*inv(ALHS)*(s1 + bc1);
    AA = I - dt*(ALHS \ A);
    b = f0 + dt*(ALHS \ (s0 + bc0));
else
    AA = I - dt*A;
    b = f0 + dt*(s0 + bc0);
end

f1 = AA\b; % Solve at each timestep

end


%% Crank Nicolson (second order implicit time advance)
function f1 = crank_nicolson(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

s0 = source_vector(pde,opts,hash_table,t);
s1 = source_vector(pde,opts,hash_table,t+dt);

bc0 = boundary_condition_vector(pde,opts,hash_table,t);
bc1 = boundary_condition_vector(pde,opts,hash_table,t+dt);

if opts.time_independent_A
    persistent A;
    persistent AA_inv;
    persistent ALHS_inv;
else
    A = [];
    AA_inv = [];
    ALHS_inv = [];
end

if isempty(AA_inv) || ~opts.time_independent_A
    
    applyLHS = ~isempty(pde.termsLHS);
        
    if applyLHS
        [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
        I = eye(numel(diag(A)));
        ALHS_inv = inv(ALHS);
%         AA = 2*I - dt*(ALHS \ A);
%         b = 2*f0 + dt*(ALHS \ A)*f0 + dt*(ALHS \ (s0+s1+bc0+bc1));
        AA = 2*I - dt*(ALHS_inv * A);
        b = 2*f0 + dt*(ALHS_inv * A)*f0 + dt*(ALHS_inv * (s0+s1+bc0+bc1));
    else
        [~,A] = apply_A(pde,opts,A_data,f0,deg);
        I = eye(numel(diag(A)));
        AA = 2*I - dt*A;
        b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
    end
    
    if ~opts.quiet; disp(['    rcond(AA) : ', num2str(rcond(AA))]); end
    
    AA_inv = inv(AA);
    
%     f1 = AA\b; % Solve at each timestep
    f1 = AA_inv * b;
else
    applyLHS = ~isempty(pde.termsLHS);
    
    if applyLHS
        I = eye(numel(diag(A)));
        AA = 2*I - dt*(ALHS_inv * A);
        b = 2*f0 + dt*(ALHS_inv * A)*f0 + dt*(ALHS_inv * (s0+s1+bc0+bc1));
    else
        I = eye(numel(diag(A)));
        AA = 2*I - dt*A;
        b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
    end
    
    f1 = AA_inv * b;
end

end
