function f = time_advance(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax)

if strcmp(opts.timestep_method,'BE')
    % Backward Euler (BE) first order
    f = backward_euler(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'FE')
    % Forward Euler (FE) first order
    f = forward_euler(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'SSPRK2')
    % Strong stability-preserving Runge-Kutta Method (also Heun's Method)
    f = ssprk2(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'SSPRK3')
    % Strong stability-preserving Runge-Kutta Method (also Heun's Method)
    f = ssprk3(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'matrix_exponential')
    % Matrix Exponential (ME) all order
    f = matrix_exponential(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'time_independent')
    % time independent d/dt==0
    f = time_independent(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif strcmp(opts.timestep_method,'CN')
    % Crank Nicolson (CN) second order
    f = crank_nicolson(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
elseif sum(strcmp(opts.timestep_method,{'ode15i','ode15s','ode45','ode23s'}))>0
    % One of the atlab ODE integrators.
    f = ODEm(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
else
    % RK3 TVD
    f = RungeKutta3(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
end

end

%% Use Matlab ODE solvers

function fval = ODEm(pde,opts,A_data,f0,t0,dt,deg,hash_table,Vmax,Emax)

clear strCR;

applyLHS = ~isempty(pde.termsLHS);

%     function res = fun_for_jacobianest(x)
%     res = explicit_ode(t0,x);
%     disp('calling ode');
%     end

    function dfdt = explicit_ode(t,f)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
        if applyLHS
            [f1,~,ALHS] = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
            dfdt = ALHS\(f1 + source + bc);       
        else
            f1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
            dfdt = f1 + source + bc;       
        end
    end

    function res = implicit_ode(t,f,dfdt)
        bc = boundary_condition_vector(pde,opts,hash_table,t);
        source = source_vector(pde,opts,hash_table,t);
        f1 = apply_A(pde,opts,A_data,f,deg,Vmax,Emax);
        res = dfdt - (f1 + source + bc);
    end

if strcmp(opts.timestep_method,'ode45')
    
    stats = 'off';
    if(~opts.quiet)
        disp('Using ode45');
        stats = 'on';
    end
    options = odeset('RelTol',1e-3,'AbsTol',1e-6,'Stats',stats);
    [tout,fout] = ode45(@explicit_ode,[t0 t0+dt],f0,options);
    
elseif strcmp(opts.timestep_method,'ode15s')
    
    %     % estimate Jacobian
    %     numjacopts.diffvar = 2;
    %     numjacopts.vectvars = [];
    %     numjacopts.thresh = 1e-10;
    %     numjacopts.fac = [];
    
    %     disp('running odenumjac')
    %     rhs = feval(@explicit_ode,t0,f0);
    %     J = odenumjac(@explicit_ode,{t0 f0},rhs,numjacopts);
    %     S = sparse(J~=0.0);
    %     disp('done')
    
    %     disp('running jacobianest')
    %     [J2,err] = jacobianest(@fun_for_jacobianest,f0);
    %     disp('done');
    
    % call ode15s
    stats = 'off';
    output_func = '';
    if(~opts.quiet)
        disp('Using ode15s');
        stats = 'on';
        output_func = @odetpbar;
    end
    options = odeset('RelTol',1e-6,'AbsTol',1e-8,...
        'Stats',stats,'OutputFcn',output_func,'Refine',20);%,'Jacobian', J2);%'JPattern',S);
    [tout,fout] = ode15s(@explicit_ode,[t0 t0+dt],f0,options);
    
elseif strcmp(opts.timestep_method,'ode23s')
    
    % call ode23s
    stats = 'off';
    output_func = '';
    if(~opts.quiet)
        disp('Using ode23s');
        stats = 'on';
        output_func = @odetpbar;
    end
    options = odeset('Stats',stats,'OutputFcn',output_func);
    [tout,fout] = ode23s(@explicit_ode,[t0 t0+dt],f0,options);
    
elseif strcmp(opts.timestep_method,'ode15i')

    if applyLHS
        error('ERROR: ode15i not yet implemented for LHS=true');
    end
    
    dfdt0 = f0.*0;
    [f0,dfdt0,resnrm] = decic(@implicit_ode,t0,f0,f0.*0+1,dfdt0,[]);
    if(~opts.quiet);disp('Using ode15i');end
    [tout,fout] = ode15i(@implicit_ode,[t0 t0+dt],f0,dfdt0);
    
end

fval = reshape(fout(end,:),size(f0));

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
%     f1 = f0 + dt * (apply_A(pde,opts,A_data,f0,deg) + s0 + bc0);
    f1 = cell(numel(f0),1);
    rhs = pde.nlinterms{1}(pde,opts,f0);
    for i = 1 : numel(f0)
        f1{i} = f0{i} + dt * rhs{i};
    end
end

end

%% SSPRK2
function f1 = ssprk2(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

if applyLHS
    error('apply LHS not implemented for SSPRK2');
else
    f1 = cell(numel(f0),1);
    rhs0 = pde.nlinterms{1}(pde,opts,f0);
    for i = 1 : numel(f0)
        f1{i} = f0{i} + dt * rhs0{i};
    end
    rhs1 = pde.nlinterms{1}(pde,opts,f1);
    for i = 1 : numel(f0)
        f1{i} = f0{i} + 0.5 * dt * ( rhs0{i} + rhs1{i} );
    end
end

end

%% SSPRK3
function f1 = ssprk3(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

if applyLHS
    error('apply LHS not implemented for SSPRK3');
else
    f1 = cell(numel(f0),1);
    rhs0 = pde.nlinterms{1}(pde,opts,f0);
    for i = 1 : numel(f0)
        f1{i} = f0{i} + dt * rhs0{i};
    end
    rhs1 = pde.nlinterms{1}(pde,opts,f1);
    for i = 1 : numel(f0)
        f1{i} = f0{i} + 0.25 * dt * ( rhs0{i} + rhs1{i} );
    end
    rhs2 = pde.nlinterms{1}(pde,opts,f1);
    for i = 1 : numel(f0)
        f1{i} = f0{i} + dt * ( rhs0{i} + rhs1{i} + 4.0 * rhs2{i} ) / 6.0;
    end
end

end

%% Backward Euler (first order implicit time advance)
function f1 = backward_euler(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

s0 = source_vector(pde,opts,hash_table,t+dt);
bc0 = boundary_condition_vector(pde,opts,hash_table,t+dt);

fun = @(g) fsolve_wrapper(pde,opts,f0,g,dt);

f0_vec = zeros(numel(f0)*numel(f0{1}),1);
for i = 1 : numel(f0)
    f0_vec((i-1)*numel(f0{1})+1:i*numel(f0{1}))=f0{i};
end

f1_vec = fsolve(fun,f0_vec);

f1 = cell(numel(f0),1);
for i = 1 : numel(f0)
    f1{i}=f1_vec((i-1)*numel(f0{1})+1:i*numel(f0{1}));
end

% if opts.time_independent_A || opts.time_independent_build_A
%     persistent A;
% end
% 
% if opts.time_independent_A % relies on inv() so is no good for poorly conditioned problems
%     persistent AA_inv;
%     persistent ALHS_inv;
%     
%     if isempty(AA_inv)
%         if applyLHS
%             [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
%             I = speye(numel(diag(A)));
%             ALHS_inv = inv(ALHS);
%             AA = I - dt*(ALHS_inv * A);
%             b = f0 + dt*(ALHS_inv * (s0 + bc0));
%         else
%             [~,A] = apply_A(pde,opts,A_data,f0,deg);
%             I = speye(numel(diag(A)));
%             AA = I - dt*A;
%             b = f0 + dt*(s0 + bc0);
%         end
%         
%         if numel(AA(:,1)) <= 4096
%             condAA = condest(AA);
%             if ~opts.quiet; disp(['    condest(AA) : ', num2str(condAA,'%.1e')]); end
%             
%             if condAA > 1e6
%                 disp(['WARNING: Using time_independent_A=true for poorly conditioned system not recommended']);
%                 disp(['WARNING: cond(A) = ', num2str(condAA,'%.1e')]);
%             end
%         end
%         
%         AA_inv = inv(AA);
%         f1 = AA_inv * b;
%     else
%         if applyLHS
%             b = f0 + dt*(ALHS_inv * (s0 + bc0));
%         else
%             b = f0 + dt*(s0 + bc0);
%         end
%         f1 = AA_inv * b;
%     end
%     
% else % use the backslash operator instead
%     if applyLHS
%         if opts.time_independent_build_A
%             if isempty(A);[~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);end
%         else
%             [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
%         end
%         I = speye(numel(diag(A)));
%         AA = I - dt*(ALHS \ A);
%         b = f0 + dt*(ALHS \ (s0 + bc0));
%     else
%         if opts.time_independent_build_A
%             if isempty(A);[~,A] = apply_A(pde,opts,A_data,f0,deg);end
%         else
%             [~,A] = apply_A(pde,opts,A_data,f0,deg);
%         end
%         I = speye(numel(diag(A)));
%         AA = I - dt*A;
%         b = f0 + dt*(s0 + bc0);
%     end
%     
%     % Direct solve
%     rescale = false;
%     if rescale
%         %     [AA_rescaled,diag_scaling] = rescale2(AA);
%         %     AA_thresholded = sparsify(AA_rescaled,1e-5);
%         [P,R,C] = equilibrate(AA);
%         AA_rescaled = R*P*AA*C;
%         b_rescaled = R*P*b;
%         f1 = AA_rescaled \ b_rescaled;
%         f1 = C*f1;
%     else
%         f1 = AA \ b;
%     end
%     
%     %     % Pre-kron solve
%     %     num_terms = numel(pde.terms);
%     %     num_dims = numel(pde.dimensions);
%     %     for d=1:num_dims
%     %         dim_mat_list{d} = zeros(size(pde.terms{1}.terms_1D{1}.mat));
%     %         for t=1:num_terms
%     %             dim_mat_list{d} = dim_mat_list{d} + pde.terms{t}.terms_1D{d}.mat;
%     %         end
%     %     end
%     %     for d=1:num_dims
%     %         A = dim_mat_list{d};
%     %         I = speye(numel(diag(A)));
%     %         AA_d =  I - dt*A;
%     %         dim_mat_inv_list{d} = inv(AA_d);
%     %     end
%     %     use_kronmultd = true;
%     %     if use_kronmultd
%     %         f1a = kron_multd(num_dims,dim_mat_inv_list,b);
%     %     else
%     %         f1a = kron_multd_full(num_dims,dim_mat_inv_list,b);
%     %     end
%     
%     %     % Iterative solve
%     %     restart = [];
%     %     tol=1e-6;
%     %     maxit=1000;
%     %     tic;
%     %     [f10,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit);
%     %     t_gmres = toc;
%     %     figure(67)
%     %     semilogy(resvec);
%     %
%     %     % Direct solve - LU approach to reducing repeated work
%     %     [L,U,P] = lu(AA);
%     %     n=size(AA,1);
%     %     ip = P*reshape(1:n,n,1); % generate permutation vector
%     %     err = norm(AA(ip,:)-L*U,1); % err should be small
%     %     tol = 1e-9;
%     %     isok = (err <= tol * norm(AA,1) ); % just a check
%     %     disp(['isok for LU: ', num2str(isok)]);
%     %     % to solve   a linear system    A * x = b, we have   P * A * x = P*b
%     %     % then from LU factorization, we have (L * U) * x = (P*b),  so    x  = U \ (L \ ( P*b))
%     %     f1_LU =   U \ (L \  (P*b));
%     %
%     %     % Direct solve - QR reduced rank
%     %     n=size(AA,1);
%     %     [Q,R,P] = qr(AA); % A*P = Q*R
%     %     % where R is upper triangular,   Q is orthogonal, Q’*Q is identity, P is column permutation
%     %     err = norm( AA*P - Q*R,1);
%     %     tol=1e-9;
%     %     isok = (err <= tol * norm(AA,1));  % just a check
%     %     disp(['isok for QR: ', num2str(isok)]);
%     %     % to solve   A * x = b,   we have A * P * (P’*x) = b, (Q*R) * (P’*x) = b
%     %     % y =   R\(Q’*b),   P*y = x
%     %     tol=1e-9;
%     %     is_illcond = abs(R(n,n)) <= tol * abs(R(1,1));
%     %     if(is_illcond)
%     %         disp('is_illcond == true');
%     %         bhat = Q'*b;
%     %         y = zeros(n,1);
%     %         yhat =  R(1:(n-1), 1:(n-1)) \ bhat(1:(n-1));  % solve with smaller system
%     %         y(1:(n-1)) = yhat(1:(n-1));
%     %         y(n) = 0;  % force last component to be zero
%     %         f1_QR = P * y;
%     %     else
%     %         disp('is_illcond == false');
%     %         f1_QR = P * (R \ (Q'*b));
%     %     end
%     %
%     %     disp(['f1-f10:  ',num2str(norm(f1-f10)/norm(f1))]);
%     %     disp(['f1-f1_LU:  ',num2str(norm(f1-f1_LU)/norm(f1))]);
%     %     disp(['f1-f1_QR:  ',num2str(norm(f1-f1_QR)/norm(f1))]);
%     %     disp(['direct runtime: ', num2str(t_direct)]);
%     %     disp(['gmres runtime: ', num2str(t_gmres)]);
%     %
%     %     f1 = f1_QR;
%     
% end
end

function [ Fg ] = fsolve_wrapper(pde,opts,f0,g,dt)

rho=cell(numel(f0),1);
for i = 1 : numel(f0)
    rho{i} = g((i-1)*numel(f0{1})+1:i*numel(f0{1}));
end
rhs = pde.nlinterms{1}(pde,opts,rho);

Fg = zeros(size(g));

for i = 1 : numel(f0)
    Fg((i-1)*numel(f0{1})+1:i*numel(f0{1}))...
        = g((i-1)*numel(f0{1})+1:i*numel(f0{1}))-dt.*rhs{i}-f0{i};
end

end

%% Crank Nicolson (second order implicit time advance)
function f1 = crank_nicolson(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

applyLHS = ~isempty(pde.termsLHS);

s0 = source_vector(pde,opts,hash_table,t);
s1 = source_vector(pde,opts,hash_table,t+dt);

bc0 = boundary_condition_vector(pde,opts,hash_table,t);
bc1 = boundary_condition_vector(pde,opts,hash_table,t+dt);

if opts.time_independent_A % uses inv() so no good for poorly conditioned systems
    persistent A;
    persistent AA_inv;
    persistent ALHS_inv;
    
    if isempty(AA_inv)
        if applyLHS
            [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            ALHS_inv = inv(ALHS);
            AA = 2*I - dt*(ALHS_inv * A);
            b = 2*f0 + dt*(ALHS_inv * A)*f0 + dt*(ALHS_inv * (s0+s1+bc0+bc1));
        else
            [~,A] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            AA = 2*I - dt*A;
            b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
        end
        
        rcondAA = rcond(AA);
        if ~opts.quiet; disp(['    rcond(AA) : ', num2str(rcondAA)]); end
        
        if 1/rcondAA > 1e6
            disp(['WARNING: Using time_independent_A=true for poorly conditioned system not recommended']);
            disp(['WARNING: cond(A) = ', num2str(rcondAA)]);
        end
        
        AA_inv = inv(AA);
        f1 = AA_inv * b;
    else
        if applyLHS
            b = 2*f0 + dt*(ALHS_inv * A)*f0 + dt*(ALHS_inv * (s0+s1+bc0+bc1));
        else
            b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
        end
        f1 = AA_inv * b;
    end
    
else % use the backslash operator for time_indepent_A = false
    if applyLHS
        [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
        I = speye(numel(diag(A)));
        AA = 2*I - dt*(ALHS \ A);
        b = 2*f0 + dt*(ALHS \ A)*f0 + dt*(ALHS \ (s0+s1+bc0+bc1));
    else
        [~,A] = apply_A(pde,opts,A_data,f0,deg);
        I = speye(numel(diag(A)));
        AA = 2*I - dt*A;
        b = 2*f0 + dt*A*f0 + dt*(s0+s1) + dt*(bc0+bc1);
    end
    f1 = AA \ b;
end

end
