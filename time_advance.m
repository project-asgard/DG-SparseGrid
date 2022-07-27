function [f,sol,hash_table] = time_advance(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax)

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
    [f,hash_table] = imex(pde,opts,A_data,f,t,dt,deg,hash_table,Vmax,Emax);
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

% when using the matlab ODE integrators we just pass them the start and
% final time (when 'time_independent_A','true') and have it return the
% intermeadiate times

if opts.time_independent_build_A
    tf = (opts.num_steps+1) * dt; % the +1 here is to account for very small deviations in the end time requested from sol
    times = [t0:dt:tf];
else
    times = [t0,t0+dt];
end

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
    sol = ode45(@explicit_ode,[t0 t0+dt],f0,options);
    fout = deval(sol,t0+dt);
    
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
%      [tout,fout0] = ode15s(@explicit_ode,[t0 t0+dt],f0,options);
     sol = ode15s(@explicit_ode,times,f0,options);
     fout = deval(sol,t0+dt);
    
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
    sol = ode23s(@explicit_ode,[t0 t0+dt],f0,options);
    fout = deval(sol,t0+dt);
    
elseif strcmp(opts.timestep_method,'ode15i')

    if applyLHS
        error('ERROR: ode15i not yet implemented for LHS=true');
    end
    
    dfdt0 = f0.*0;
    [f0,dfdt0,resnrm] = decic(@implicit_ode,t0,f0,f0.*0+1,dfdt0,[]);
    if(~opts.quiet);disp('Using ode15i');end
    sol = ode15i(@implicit_ode,[t0 t0+dt],f0,dfdt0);
    fout = deval(sol,t0+dt);
end

% fval = reshape(fout(end,:),size(f0));

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

if opts.time_independent_A || opts.time_independent_build_A
    persistent A;
end

if opts.time_independent_A % relies on inv() so is no good for poorly conditioned problems
    persistent AA_inv;
    persistent ALHS_inv;
    
    if isempty(AA_inv)
        if applyLHS
            [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            ALHS_inv = inv(ALHS);
            AA = I - dt*(ALHS_inv * A);
            b = f0 + dt*(ALHS_inv * (s0 + bc0));
        else
            [~,A] = apply_A(pde,opts,A_data,f0,deg);
            I = speye(numel(diag(A)));
            AA = I - dt*A;
            b = f0 + dt*(s0 + bc0);
        end
        
        if numel(AA(:,1)) <= 4096
            condAA = condest(AA);
            if ~opts.quiet; disp(['    condest(AA) : ', num2str(condAA,'%.1e')]); end
            
            if condAA > 1e6
                disp(['WARNING: Using time_independent_A=true for poorly conditioned system not recommended']);
                disp(['WARNING: cond(A) = ', num2str(condAA,'%.1e')]);
            end
        end
        
        AA_inv = inv(AA);
        f1 = AA_inv * b;
    else
        if applyLHS
            b = f0 + dt*(ALHS_inv * (s0 + bc0));
        else
            b = f0 + dt*(s0 + bc0);
        end
        f1 = AA_inv * b;
    end
    
else % use the backslash operator instead
    if applyLHS
        if opts.time_independent_build_A
            if isempty(A);[~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);end
        else
            [~,A,ALHS] = apply_A(pde,opts,A_data,f0,deg);
        end
        I = speye(numel(diag(A)));
        AA = I - dt*(ALHS \ A);
        b = f0 + dt*(ALHS \ (s0 + bc0));
    else
        if opts.time_independent_build_A
            if isempty(A);[~,A] = apply_A(pde,opts,A_data,f0,deg);end
        else
            [~,A] = apply_A(pde,opts,A_data,f0,deg);
        end
        I = speye(numel(diag(A)));
        AA = I - dt*A;
        b = f0 + dt*(s0 + bc0);
    end
    
    % Direct solve
    rescale = false;
    if rescale
        %     [AA_rescaled,diag_scaling] = rescale2(AA);
        %     AA_thresholded = sparsify(AA_rescaled,1e-5);
        [P,R,C] = equilibrate(AA);
        AA_rescaled = R*P*AA*C;
        b_rescaled = R*P*b;
        f1 = AA_rescaled \ b_rescaled;
        f1 = C*f1;
    else
        f1 = AA \ b;
    end
    
    %     % Pre-kron solve
    %     num_terms = numel(pde.terms);
    %     num_dims = numel(pde.dimensions);
    %     for d=1:num_dims
    %         dim_mat_list{d} = zeros(size(pde.terms{1}.terms_1D{1}.mat));
    %         for t=1:num_terms
    %             dim_mat_list{d} = dim_mat_list{d} + pde.terms{t}.terms_1D{d}.mat;
    %         end
    %     end
    %     for d=1:num_dims
    %         A = dim_mat_list{d};
    %         I = speye(numel(diag(A)));
    %         AA_d =  I - dt*A;
    %         dim_mat_inv_list{d} = inv(AA_d);
    %     end
    %     use_kronmultd = true;
    %     if use_kronmultd
    %         f1a = kron_multd(num_dims,dim_mat_inv_list,b);
    %     else
    %         f1a = kron_multd_full(num_dims,dim_mat_inv_list,b);
    %     end
    
    %     % Iterative solve
    %     restart = [];
    %     tol=1e-6;
    %     maxit=1000;
    %     tic;
    %     [f10,flag,relres,iter,resvec] = gmres(AA,b,restart,tol,maxit);
    %     t_gmres = toc;
    %     figure(67)
    %     semilogy(resvec);
    %
    %     % Direct solve - LU approach to reducing repeated work
    %     [L,U,P] = lu(AA);
    %     n=size(AA,1);
    %     ip = P*reshape(1:n,n,1); % generate permutation vector
    %     err = norm(AA(ip,:)-L*U,1); % err should be small
    %     tol = 1e-9;
    %     isok = (err <= tol * norm(AA,1) ); % just a check
    %     disp(['isok for LU: ', num2str(isok)]);
    %     % to solve   a linear system    A * x = b, we have   P * A * x = P*b
    %     % then from LU factorization, we have (L * U) * x = (P*b),  so    x  = U \ (L \ ( P*b))
    %     f1_LU =   U \ (L \  (P*b));
    %
    %     % Direct solve - QR reduced rank
    %     n=size(AA,1);
    %     [Q,R,P] = qr(AA); % A*P = Q*R
    %     % where R is upper triangular,   Q is orthogonal, Q’*Q is identity, P is column permutation
    %     err = norm( AA*P - Q*R,1);
    %     tol=1e-9;
    %     isok = (err <= tol * norm(AA,1));  % just a check
    %     disp(['isok for QR: ', num2str(isok)]);
    %     % to solve   A * x = b,   we have A * P * (P’*x) = b, (Q*R) * (P’*x) = b
    %     % y =   R\(Q’*b),   P*y = x
    %     tol=1e-9;
    %     is_illcond = abs(R(n,n)) <= tol * abs(R(1,1));
    %     if(is_illcond)
    %         disp('is_illcond == true');
    %         bhat = Q'*b;
    %         y = zeros(n,1);
    %         yhat =  R(1:(n-1), 1:(n-1)) \ bhat(1:(n-1));  % solve with smaller system
    %         y(1:(n-1)) = yhat(1:(n-1));
    %         y(n) = 0;  % force last component to be zero
    %         f1_QR = P * y;
    %     else
    %         disp('is_illcond == false');
    %         f1_QR = P * (R \ (Q'*b));
    %     end
    %
    %     disp(['f1-f10:  ',num2str(norm(f1-f10)/norm(f1))]);
    %     disp(['f1-f1_LU:  ',num2str(norm(f1-f1_LU)/norm(f1))]);
    %     disp(['f1-f1_QR:  ',num2str(norm(f1-f1_QR)/norm(f1))]);
    %     disp(['direct runtime: ', num2str(t_direct)]);
    %     disp(['gmres runtime: ', num2str(t_gmres)]);
    %
    %     f1 = f1_QR;
    
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

function [f1,hash_table] = imex(pde,opts,A_data,f0,t,dt,deg,hash_table,Vmax,Emax)

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

persistent hash_table_FG A_data_FG perm_FG iperm_FG pvec_FG
persistent B B0
persistent hash 
persistent x v FMWT_x
persistent output
%persistent perm_SG iperm_SG pvec_SG
%persistent per iper %Conversion for limiters
%persistent FMWT_2D
%persistent B0 B M FMWT_x x

%Switch for IMEX iteration: 
% 1 : Backward Euler in explicit terms 'E',
%     Backward Euler in implicit terms 'I'.
% 0 : SSP-RK2 IMEX.  'E' terms are actually handled explicitly.

%%% C++ : Focus on the case BEFE = 0

BEFE = 0;

pos_adapt = true;
adapt = true;
limit = true;
pos_tol = -1e-14;
lim_tol =  1e-14;

%Get quadrature points in realspace stiffness matrix calculation
[Meval,nodes] = matrix_plot_D(pde,opts,pde.dimensions{1});

%Make fake pde file for use in physical realspace transformation
pde_1d.dimensions = pde.dimensions(1);

if isempty(hash_table_FG)
    %hash_table_FG = hash_table_nD(pde.get_lev_vec(),'FG');
    [hash_table_FG.elements, hash_table_FG.elements_idx]  = hash_table_sparse_nD (pde.get_lev_vec, opts.max_lev, 'FG');
    A_data_FG = global_matrix(pde,opts,hash_table_FG);
    [perm_FG,iperm_FG,pvec_FG] = sg_to_fg_mapping_2d(pde,opts,A_data_FG);
    x = pde.dimensions{1}.min:(pde.dimensions{1}.max-pde.dimensions{1}.min)/2^pde.dimensions{1}.lev:pde.dimensions{1}.max;
    v = pde.dimensions{2}.min:(pde.dimensions{2}.max-pde.dimensions{2}.min)/2^pde.dimensions{2}.lev:pde.dimensions{2}.max;
    FMWT_x = OperatorTwoScale_wavelet2(opts.deg,pde.dimensions{1}.lev);    
end

moment_flag = (numel(pde.moments) == 3);



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
    
    if isempty(hash)
        hash = hash_table;
    end
    A_data = global_matrix(pde,opts,hash);
    if moment_flag
        moment_mat = cell(numel(pde.moments),1);
        for i=1:numel(pde.moments)
            moment_mat{i} = moment_reduced_matrix(opts,pde,A_data,hash_table,i);
        end
    end
    
    if isempty(B)
        [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
        
        %Positive Perserving Initial Condition
        Q = B0*f0;
        fprintf('POS:     min(B0*f0) = %e.  Negative on %d elements.\n',min(Q),sum(Q < pos_tol));
        
        if pos_adapt && 0
            pos_bool = true;
            while pos_bool
            %%% Adaptivity stuff       
                %figure(16); imagesc(QQ); title(sprintf('Negative Elements = %4d',sum(QQ(:))));
                if min(Q) < pos_tol
                    %[hash_pos,A_pos] = addNegativeElements(pde,opts,hash,Q,pos_tol);
                    [hash_pos,A_pos] = hierarchicalPostivity(pde,opts,f0,hash,A_data,pos_tol);
                    fprintf('POS:     --> Added %d dof.  ',numel(hash_pos.elements_idx)-numel(hash.elements_idx));
                else
                    pos_bool = false;
                    hash_pos = hash;
                    break
                end
                %Check to see if hash_table has changed
                if numel(hash.elements_idx) ~= numel(hash_pos.elements_idx)
                    hash = hash_pos;
                    A_data = A_pos;
                    f0 = initial_condition_vector(pde, opts, hash, t);
                    [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);

                    %Update moment matrices
                    if moment_flag
                        for i=1:numel(pde.moments)
                            moment_mat{i} = moment_reduced_matrix(opts,pde,A_data,hash,i);
                        end
                    end
                    hash_table_1D = hash_table_2D_to_1D(hash,opts);
                    f_2s = f0 + dt*fast_2d_matrix_apply(opts,pde,A_data,f0,'E');
                    Q = B0*f_2s;
                    fprintf('Updated f0: min(B0*f0) = %e.\n',min(Q));
                else
                    fprintf('Updated f0: min(B0*f0) = %e.\n',min(Q));
                    pos_bool = false;
                end

            end
        end
    end
    
    f0_hash = hash;
    figure(20); plotVec(x,v,1,B*f0,[],0);
    
    %[hash_new,A_new,~] = addLimitElements(pde,opts,hash,B*f0,f0,0.1);
    
    %Explicit step
    f_2s = f0 + dt*fast_2d_matrix_apply(opts,pde,A_data,f0,'E');
    
    %Adapt solution
    if adapt 
        [pde,hash_new,~,A_new,ele_added] = adapt_stripped(pde,opts,hash,f_2s,'r');
        if numel(ele_added) > 0
            fprintf('ADAPT:   Adaptivity on f_2s added %d elements.\n',numel(ele_added));
            f0 = fill_from_hash_table(pde,opts,hash,hash_new,f0);
            hash = hash_new;
            A_data = A_new;
            f_2s = f0 + dt*fast_2d_matrix_apply(opts,pde,A_data,f0,'E');
            [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
        end
    end
    
    if isempty(output)
        output.data_vec = numel(hash.elements_idx);
        output.min_vec = min(B0*f0);
        output.max_vec = max(B0*f0);
        output.T = 0;
    end
    
    %one = new_md_func(2,{@(x,p,t,dat) 0*x+1,@(v,t,p,dat) 0*v+1, @(t,p) 0*t+1});
    %f_2s = md_eval_function(opts, opts.deg, pde.dimensions, ...
    %pde.params, {one}, hash_table, pde.transform_blocks, t);
    dof = numel(f_2s);
    Q = B0*f_2s;
    fprintf('POS:     min(B0*f_2s) = %e.  Negative on %d elements.\n',min(Q),sum(Q < pos_tol))
    if pos_adapt
        pos_bool = true;
        while pos_bool
        %%% Adaptivity stuff       
            %QQ = B*f_2s; QQ = reshape(QQ,2^pde.dimensions{1}.lev,[]); QQ = flipud(QQ); QQ = (QQ < pos_tol);
            %figure(16); imagesc(QQ); title(sprintf('Negative Elements = %4d',sum(QQ(:))));
            if min(Q) < pos_tol
                %[hash_pos,A_pos] = addNegativeElements(pde,opts,hash,Q,pos_tol);
                [hash_pos,A_pos] = hierarchicalPostivity(pde,opts,f_2s,hash,A_data,pos_tol);
                fprintf('POS:     --> Added %d dof.  ',numel(hash_pos.elements_idx)-numel(hash.elements_idx));
            else
                pos_bool = false;
                hash_pos = hash;
                break
            end
            %Check to see if hash_table has changed
            if numel(hash.elements_idx) ~= numel(hash_pos.elements_idx)
                f0 = fill_from_hash_table(pde,opts,hash,hash_pos,f0);
                f0_hash = hash_pos;
                hash = hash_pos;
                A_data = A_pos;
                [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);

                %Update moment matrices
                if moment_flag
                    for i=1:numel(pde.moments)
                        moment_mat{i} = moment_reduced_matrix(opts,pde,A_data,hash,i);
                    end
                end
                hash_table_1D = hash_table_2D_to_1D(hash,opts);
                f_2s = f0 + dt*fast_2d_matrix_apply(opts,pde,A_data,f0,'E');
                Q = B0*f_2s;
                fprintf('Updated f_2s: min(B0*f_2s) = %e.\n',min(Q));
            else
                fprintf('Updated f_2s: min(B0*f_2s) = %e.\n',min(Q));
                pos_bool = false;
            end
            
        end
    end
    fprintf('POS:     Positivity preserving req for f_2s added a total of %d elements\n',(numel(f_2s)-dof)/opts.deg^2);
    fprintf('LIM:     Pre  limited f_2s has cell average min of %e\n',min(Q));
    
    dof = numel(f_2s);Qe = B*f_2s;
    if limit
        %while true
        Q = B*f_2s;
        [hash_new,A_new,f_2s] = addLimitElements(pde,opts,hash,Q,f_2s,1,0,lim_tol);
        %[hash_new,A_new,f_2sl] = addHierLimitElements(pde,opts,hash,A_data,f_2s,1,0);
        %fprintf('LIM:     New Elements Idx = [%d,%d]. Norm = %e\n',numel(f_2s)+1,numel(f_2sl),norm(f_2sl(numel(f_2s)+1:end)));f_2s=f_2sl;
        %fprintf('LIM:     New Elements Idx = [%d,%d]. Norm = %e\n',dof+1,numel(f_2sl),norm(f_2sl(dof+1:end)));f_2s=f_2sl;
        fprintf('LIM:     Limiter Added %d elements\n',numel(hash_new.elements_idx)-numel(hash.elements_idx));
        f0 = fill_from_hash_table(pde,opts,hash,hash_new,f0);
        hash = hash_new; A_data = A_new;
        [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
        figure(20); plotVec(x,v,1,B*f_2s,[],0);%view([0 90]);colorbar;
        %end
    end
    
    Q = B0*f_2s;
    fprintf('LIM:     Post limited f_2s has cell average min of %e\n',min(Q));
    output.data_vec = [output.data_vec;numel(hash.elements_idx)];
    output.min_vec = [output.min_vec;min(Q)];
    output.max_vec = [output.max_vec;max(Q)];
    output.T = [output.T;t+0.5*dt];
        
    %%%%%% Testing slope limiters for SG functions.  
    %%%%%%  C++ do not implement %%%%%%%%%
    %f_2sl = limitWrapper(pde,opts,perm_FG,iperm_FG,pvec_FG,perm_SG,iperm_SG,pvec_SG,per,iper,FMWT_2D,f_2s);
    %%%%%% END %%%%%%
    
    
    if moment_flag
        %Create rho_2s
        mom0 = moment_mat{1}*f_2s; %integral of (f,1)_v
        %mom0 = FMWT_x*slopeLimiter(x,deg-1,FMWT_x'*mom0,1,1);
        mom0_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom0,hash_table_1D);
        pde.params.n  = @(x) interp1(nodes,mom0_real,x,'nearest','extrap');

        mom1 = moment_mat{2}*f_2s; %integral of (f,v)_v
        %mom1 = FMWT_x*slopeLimiter(x,deg-1,FMWT_x'*mom1,1,1);
        mom1_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom1,hash_table_1D);
        pde.params.u  = @(x) interp1(nodes,mom1_real,x,'nearest','extrap')./pde.params.n(x);

        mom2 = moment_mat{3}*f_2s; %integral of (f,v^2)_v
        %mom2 = FMWT_x*slopeLimiter(x,deg-1,FMWT_x'*mom2,1,1);
        mom2_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom2,hash_table_1D);
        pde.params.th = @(x) interp1(nodes,mom2_real,x,'nearest','extrap')./pde.params.n(x) - pde.params.u(x).^2;
    end
        
            %Plot moments
%     if ~opts.quiet || 1
%         
%         fig1 = figure(1000);
%         fig1.Units = 'Normalized';
%         fig1.Position = [0.5 0.5 0.3 0.3];
%         subplot(2,2,1);
%         plot(nodes,pde.params.n(nodes));
%         title('n_f');
%         subplot(2,2,2);
%         plot(nodes,pde.params.u(nodes));
%         title('u_f');
%         subplot(2,2,3);
%         plot(nodes,pde.params.th(nodes));
%         title('th_f');
%         sgtitle("Fluid Variables. t = "+num2str(t+dt));
%         drawnow
%         
%     end
        
    %Update coefficients
    pde = get_coeff_mats(pde,opts,t,0);
    
    %f2 now
    [f_2,flag,relres,iter1] = bicgstabl(@(x) x - dt*fast_2d_matrix_apply(opts,pde,A_data,x,'I'),f_2s,1e-12,numel(f0));
    %[f_2,flag,relres,iter1] = gmres(@(x) x - dt*fast_2d_matrix_apply(opts,pde,A_data,x,'I'),f_2s,10,1e-12,numel(f_2s)/10);
    if flag ~= 0
        fprintf('BICGSTABL did not converge.  flag = %d, relres = %5.4e\n',flag,relres);
        assert(relres < 1e-10)
    end
    
    %%%%%
    %%% Third stage
    %%%%%
    
    f_3s = f0 + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f0+f_2,'E') ...
              + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f_2,'I');
    
    %Adapt solution
    if adapt 
        [pde,hash_new,~,A_new,ele_added] = adapt_stripped(pde,opts,hash,f_3s,'r');
        if numel(ele_added) > 0
            fprintf('ADAPT:   Adaptivity on f_3s added %d elements.\n',numel(ele_added));
            f0 = fill_from_hash_table(pde,opts,hash,hash_new,f0);
            f_2 = fill_from_hash_table(pde,opts,hash,hash_new,f_2);
            hash = hash_new;
            A_data = A_new;
            f_3s = f0 + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f0+f_2,'E') ...
                      + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f_2,'I');
            [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
        end
    end
    
    dof = numel(f_3s);
    Q = B0*f_3s;
    fprintf('POS:     min(B0*f_3s) = %e. Negative on %d elements.\n',min(Q),sum(Q < pos_tol));
    if pos_adapt
        pos_bool = true;
        while pos_bool
            if min(Q) < pos_tol
                %[hash_pos,A_pos] = addNegativeElements(pde,opts,hash,Q,pos_tol);
                [hash_pos,A_pos] = hierarchicalPostivity(pde,opts,f_3s,hash,A_data,pos_tol);
                fprintf('POS:     ---> Added %d elements.  ',numel(hash_pos.elements_idx)-numel(hash.elements_idx));
            else
                hash_pos = hash;
                pos_bool = false;
                break
            end
            %Check to see if hash_table has changed
            if numel(hash.elements_idx) ~= numel(hash_pos.elements_idx)
                f_2 = fill_from_hash_table(pde,opts,hash,hash_pos,f_2);
                f0  = fill_from_hash_table(pde,opts,hash,hash_pos,f0);
                hash = hash_pos;
                A_data = A_pos;
                [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);

                %Update moment matrices
                if moment_flag
                    for i=1:numel(pde.moments)
                        moment_mat{i} = moment_reduced_matrix(opts,pde,A_data,hash,i);
                    end
                end
                hash_table_1D = hash_table_2D_to_1D(hash,opts);
                f_3s = f0  + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f0+f_2,'E') ...
                           + 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,f_2,'I');
                Q = B0*f_3s;
                fprintf('Updated f_3s: min(B0*f_3s) = %e.\n',min(Q));
            else
                fprintf('Updated f_3s: min(B0*f_3s) = %e.\n',min(Q));
                pos_bool = false;
            end
        end
    end
    fprintf('Positivity preserving req for f_3s added a total of %d elements\n',(numel(f_3s)-dof)/opts.deg^2);
    fprintf('LIM:      Pre  limited f_3s has cell average min of %e\n',min(Q));
    
    if limit
        Q = B*f_3s;
        [hash_new,A_new,f_3s] = addLimitElements(pde,opts,hash,Q,f_3s,1,0,lim_tol);
        %fprintf('LIM:     New Elements Norm = %e\n',norm(f_3sl(numel(f_3s)+1:end)));f_3s=f_3sl;
        fprintf('LIM:     Limiter Added %d elements\n',numel(hash_new.elements_idx)-numel(hash.elements_idx));
        hash = hash_new; A_data = A_new;
        [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
        figure(20); plotVec(x,v,1,B*f_3s,[],0);
    end
    
    Q = B0*f_3s;
    fprintf('LIM:     Post limited f_3s has cell average min of %e\n',min(Q));
    output.data_vec = [output.data_vec;numel(hash.elements_idx)];
    output.min_vec = [output.min_vec;min(Q)];
    output.max_vec = [output.max_vec;max(Q)];
    output.T = [output.T;t+dt];
    
    if moment_flag
        %Create rho_3s
        mom0 = moment_mat{1}*f_3s; %integral of (f,1)_v
        %mom0 = FMWT_x*slopeLimiter(x,deg-1,FMWT_x'*mom0,1,1);
        mom0_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom0,hash_table_1D);
        pde.params.n  = @(x) interp1(nodes,mom0_real,x,'nearest','extrap');

        mom1 = moment_mat{2}*f_3s; %integral of (f,v)_v
        %mom1 = FMWT_x*slopeLimiter(x,deg-1,FMWT_x'*mom1,1,1);
        mom1_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom1,hash_table_1D);
        pde.params.u  = @(x) interp1(nodes,mom1_real,x,'nearest','extrap')./pde.params.n(x);

        mom2 = moment_mat{3}*f_3s; %integral of (f,v^2)_v
        %mom2 = FMWT_x*slopeLimiter(x,deg-1,FMWT_x'*mom2,1,1);
        mom2_real = wavelet_to_realspace(pde_1d,opts,{Meval},mom2,hash_table_1D);
        pde.params.th = @(x) interp1(nodes,mom2_real,x,'nearest','extrap')./pde.params.n(x) - pde.params.u(x).^2;
    end
            
    %Update coefficients
    pde = get_coeff_mats(pde,opts,t,0);

    [f_3,flag,relres,iter2] = bicgstabl(@(x) x - 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,x,'I'),f_3s,1e-12,numel(f_3s));
    %[f_3,flag,relres,iter2] = gmres(@(x) x - 0.5*dt*fast_2d_matrix_apply(opts,pde,A_data,x,'I'),f_3s,10,1e-12,numel(f_3s)/10);
    if flag ~= 0
        fprintf('BICGSTABL did not converge.  flag = %d, relres = %5.4e\n',flag,relres);
        assert(relres < 1e-10)
    end
    
    %Update timestep to final stage
    f1 = f_3;
    
    if moment_flag
        time = t+dt;
        if opts.case_ == 2 && (abs(pde.params.nu) < 1e-8)
            %Get analytic solution at t_{n+1}

            fval_analytic = exact_solution_vector(pde,opts,hash_table,time);
            %fval_analytic = exact_solution_vector(pde,opts,hash_table,t);

            %Get analytic moments

            %n 
            mom0_a = moment_mat{1}*fval_analytic;
            mom0_a_r = wavelet_to_realspace(pde_1d,opts,{Meval},mom0_a,hash_table_1D);
            analytic_moments.n = mom0_a_r;

            %u
            mom1_a = moment_mat{2}*fval_analytic;
            mom1_a_r = wavelet_to_realspace(pde_1d,opts,{Meval},mom1_a,hash_table_1D);
            analytic_moments.u = mom1_a_r./mom0_a_r;

            %T
            mom2_a = moment_mat{3}*fval_analytic;
            mom2_a_r = wavelet_to_realspace(pde_1d,opts,{Meval},mom2_a,hash_table_1D);
            analytic_moments.T = mom2_a_r./mom0_a_r - (analytic_moments.u).^2;

            fig1 = figure(1000);
            fig1.Units = 'Normalized';
            fig1.Position = [0.5 0.5 0.3 0.3];
            subplot(2,2,1);
            plot(nodes,pde.params.n(nodes)-analytic_moments.n');
            title('n_f');
            subplot(2,2,2);
            plot(nodes,pde.params.u(nodes)-analytic_moments.u');
            title('u_f');
            subplot(2,2,3);
            plot(nodes,pde.params.th(nodes)-analytic_moments.T');
            title('th_f');
            sgtitle("Fluid Variables. t = "+num2str(time));
            drawnow
        else
            %Plot moments

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
            sgtitle("Fluid Variables. t = "+num2str(time));
            drawnow

        end
    end
    
    if adapt
        fprintf('         Pre-coarsened element count: %d\n',numel(hash.elements_idx));
        [pde,hash_new,f1_cull,A_data,elements_dropped] = adapt_stripped(pde,opts,hash,f1,'c');
        %Check to see if hash has changed
        if numel(hash_new.elements_idx) ~= numel(hash.elements_idx)
            fprintf('Coarsened mesh by %d elements\n',numel(hash.elements_idx)-numel(hash_new.elements_idx));
            %hash_bc = hash;
            f1_bc = f1;
            hash = hash_new;
            f1 = f1_cull;
            [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
        end
        Q = B0*f1;
        num_ele = numel(hash.elements_idx);
        fprintf('POS:     min(B0*f1) = %e. Negative on %d elements.\n',min(Q),sum(Q < pos_tol));
        if pos_adapt
            pos_bool = (min(Q) < pos_tol);
            while pos_bool
                %[hash_pos,A_pos] = addNegativeElements(pde,opts,hash,Q,pos_tol);
                [hash_pos,A_pos] = hierarchicalPostivity(pde,opts,f1,hash,A_data,pos_tol);
                [f1,hash_pos,A_pos] = fill_from_hash_table(pde,opts,hash,hash_pos,f1,f1_bc,elements_dropped,A_pos);
                fprintf('POS:     ---> Added %d elements.  ',numel(hash_pos.elements_idx)-numel(hash.elements_idx));
                if numel(hash.elements_idx) == numel(hash_pos.elements_idx)
                    pos_bool = false;
                end
                hash = hash_pos;
                A_data = A_pos;
                [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
                Q = B0*f1;
                fprintf('Updated f1: min(B0*f1) = %e.\n',min(Q));
                if min(Q) >= pos_tol
                    pos_bool = false;
                end
            end
            fprintf('POS:     min(B0*f1) = %e. Added a total of %d elements.\n',min(Q),numel(hash.elements_idx)-num_ele);
        end
    end
    
    if limit
        Q = B*f1;
        [hash_new,A_new,f1] = addLimitElements(pde,opts,hash,Q,f1,1,0,lim_tol);
        %[hash_new,A_new,f1] = addLimitElements(pde,opts,hash,Q,f1,1,0,lim_tol);
        %fprintf('LIM:     New Elements Norm = %e\n',norm(f1_l(numel(f1)+1:end)));f1=f1_l;
        fprintf('LIM:     Limiter Added %d elements\n',numel(hash_new.elements_idx)-numel(hash.elements_idx));
        hash = hash_new; A_data = A_new;
        [B,B0] = WaveletToRealspaceTransMatrix(pde,opts,A_data);
        figure(20); plotVec(x,v,1,B*f1,[],0);
    end
    
    hash_table = hash;
end

if t+dt > dt*opts.num_steps - dt/10
    []; 
end

end
