classdef OPTS
    
    properties
        
        dt_set_at_runtime = true;
        CFL = 0.01;
        dt = 0;
        quiet = false;
        grid_type = 'SG';        
        timestep_method = 'RK3';
        build_A = false;
        adapt = false;
        use_oldhash = false;
        use_oldcoeffmat = false;
        time_independent_A = false;
        many_solution_capable = false;
        max_lev = 8;
        adapt_threshold = 1e-3;
        refinement_method = 1;
        adapt_initial_condition = false;
        uniform_output = false;
        save_output = false;
        output_filename_id = '';

    end
    
end