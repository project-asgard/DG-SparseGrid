classdef OPTS
    
    properties
        
        dt_set_at_runtime = true;
        CFL = 0.01;
        dt = 0;
        quiet = false;
        grid_type = 'SG'; % 'SG', 'FG'        
        timestep_method = 'RK3'; % 'RK3', 'FE', 'BE', 'ode15s', 'ode15i', 'ode45'
        build_A = false;
        adapt = false;
        use_oldhash = false;
        use_oldcoeffmat = false;
        time_independent_A = false;
        time_independent_build_A = false;
        many_solution_capable = false;
        max_lev = 8;
        adapt_threshold = 1e-3;
        refinement_method = 1;
        adapt_initial_condition = false;
        output_grid = 'quadrature'; % 'quadrature', 'fixed', 'uniform'
        save_output = false;
        output_filename_id = '';
        plot_freq = 1;
        save_freq = 1;

    end
    
end