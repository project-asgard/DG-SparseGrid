function [ err, opts ] = asgard_run_pde_system( opts, pde_system )

    %% Set time step (may need to move inside time loop)
    dt = pde_system.set_dt( pde_system, opts.CFL );
    if opts.dt_set_at_runtime
        dt = opts.dt;
    else
        opts.dt = dt;
    end
    if ~opts.quiet; fprintf('dt = %g\n', dt ); end
    
    t = 0.0;
    
    %% Set initial conditions
    pde_system.set_initial_conditions( t );
    
    %% Write initial conditions here
    
    %% Time Loop
    count = 0;
    if ~opts.quiet; disp('Advancing time ...'); end
    for L = 1 : opts.num_steps
        
        if ~opts.quiet
            fprintf( 'Step %i of %i (t= %f)\n', L, opts.num_steps, t );
        end
        
        if ~opts.quiet; disp(['    Time step (dt=',num2str(dt),')']); end
        
        time_stepper( pde_system, t, dt );
        
        t = t + dt;
        
        count = count + 1;
        
    end
    
    err = pde_system.compute_error( t );
    
    if ~opts.quiet
        for i = 1 : numel(pde_system.unknowns)
            disp(['    unknown : ', num2str(i)]);
            disp(['    num_dof : ', num2str(pde_system.unknowns{i}.size())]);
            disp(['    wavelet space absolute err : ', num2str(err(i))]);
        end
    end
    
end