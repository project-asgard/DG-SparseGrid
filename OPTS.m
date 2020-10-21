classdef OPTS
    
    properties
        
        lev = 3; % can be set as 'lev',3 or 'lev',[2,3]
        deg = 2;
        dt_set_at_runtime = true;
        CFL = 0.01;
        dt = 0;
        quiet = false;
        grid_type = 'SG'; % 'SG', 'FG'
        timestep_method = 'RK3'; % 'RK3', 'FE', 'BE', 'ode15s', 'ode15i', 'ode45', 'time_independent'
        build_A = false;
        adapt = false;
        use_connectivity = false;
        use_oldhash = false;
        use_oldcoeffmat = false;
        use_sparse_A = false;
        time_independent_A = false;
        time_independent_build_A = false;
        many_solution_capable = false; % everything needs to transition to this by defalt, all it means is that a soln or init_cond can be defined as the sum of terms, i.e., soln1 + soln2 defined as solutions = {soln1, soln2}, where each soln is a multi-d function definition
        max_lev = 8;
        max_lev_coeffs = true; % when enabled, build partial term coeff matrices for some max level,
        % and rechain 1d term matrices as adaptivity dictates
        adapt_threshold = 1e-3;
        refinement_method = 1;
        adapt_initial_condition = false;
        output_grid = 'quadrature'; % 'quadrature', 'fixed', 'uniform', 'quadrature_with_end_points'
        save_output = false;
        output_filename_id = '';
        plot_freq = 1;
        save_freq = 1;
        case_ = 1;
        calculate_mass = true; % calculate and print the mass
        normalize_by_mass = false; % normalize the analytic solution by the initial cond mass
        start_time = 0;
        num_steps = 5;   
        
    end
    
    methods
        
        function opts = OPTS(varargin)
            
            if nargin == 0
            else
                input_parser = inputParser();
                input_parser.KeepUnmatched = true;
                
                valid_grid_types = {'SG','FG'};
                check_grid_type = @(x) any(validatestring(x,valid_grid_types));
                valid_timestep_methods = {'BE','CN','ode15i','ode15s','ode45','RK3','FE','time_independent'};
                check_timestep_method = @(x) any(strcmp(x,valid_timestep_methods));
                valid_output_grids = {'quadrature','fixed','uniform','quadrature_with_end_points'};
                check_output_grid = @(x) any(strcmp(x,valid_output_grids));
                
                addOptional(input_parser,'lev',opts.lev, @isnumeric);
                addOptional(input_parser,'deg',opts.deg, @isnumeric);
                addOptional(input_parser,'start_time',opts.start_time, @isnumeric);
                addOptional(input_parser,'num_steps',opts.num_steps, @isnumeric);
                addOptional(input_parser,'quiet',opts.quiet,@islogical);
                addOptional(input_parser,'timestep_method',opts.timestep_method, check_timestep_method);
                addOptional(input_parser,'grid_type',opts.grid_type, check_grid_type);
                addOptional(input_parser,'CFL',opts.CFL, @isnumeric);
                addOptional(input_parser,'dt',opts.dt, @isnumeric);
                addOptional(input_parser,'adapt',opts.adapt, @islogical);
                addOptional(input_parser,'use_oldhash',opts.use_oldhash, @islogical);
                addOptional(input_parser,'use_oldcoeffmat',opts.use_oldcoeffmat, @islogical);
                addOptional(input_parser,'time_independent_A',opts.time_independent_A,@islogical);
                addOptional(input_parser,'time_independent_build_A',opts.time_independent_build_A,@islogical);
                addOptional(input_parser,'many_solution_capable',opts.many_solution_capable,@islogical);
                addOptional(input_parser,'max_lev',opts.max_lev, @isnumeric);
                addOptional(input_parser,'max_lev_coeffs',opts.max_lev_coeffs, @islogical);
                addOptional(input_parser,'adapt_threshold',opts.adapt_threshold, @isnumeric);
                addOptional(input_parser,'refinement_method',opts.refinement_method, @isnumeric);
                addOptional(input_parser,'adapt_initial_condition',opts.adapt_initial_condition,@islogical);
                addOptional(input_parser,'save_output',opts.save_output,@islogical);
                addOptional(input_parser,'output_filename_id',opts.output_filename_id,@ischar);
                addOptional(input_parser,'plot_freq',opts.plot_freq, @isnumeric);
                addOptional(input_parser,'save_freq',opts.save_freq, @isnumeric);
                addOptional(input_parser,'output_grid',opts.output_grid,check_output_grid);
                addOptional(input_parser,'use_connectivity',opts.use_connectivity,@islogical);
                addOptional(input_parser,'use_sparse_A',opts.use_sparse_A,@islogical);
                addOptional(input_parser,'case',opts.case_,@isnumeric);
                addOptional(input_parser,'calculate_mass',opts.calculate_mass,@islogical);
                addOptional(input_parser,'normalize_by_mass',opts.normalize_by_mass,@islogical);
                
                parse(input_parser,varargin{:})
                
                % CFL priority
                % low        : opts.CFL
                % med        : pde.CFL
                % high       : command line CFL
                % extra high : set dt at command line
                
                CFL_set = true;
                opts.dt_set_at_runtime = true;
                for c=input_parser.UsingDefaults
                    if strcmp(c{1},'CFL')
                        CFL_set = false;
                    end
                    if strcmp(c{1},'dt')
                        opts.dt_set_at_runtime = false;
                    end
                end
                
                opts.CFL = opts.CFL;
                
                if CFL_set
                    opts.CFL = input_parser.Results.CFL;
                end
                
                opts.lev = input_parser.Results.lev;
                opts.deg = input_parser.Results.deg;
                opts.dt = input_parser.Results.dt;
                opts.quiet = input_parser.Results.quiet;
                opts.grid_type = input_parser.Results.grid_type;
                opts.timestep_method = input_parser.Results.timestep_method;
                opts.adapt = input_parser.Results.adapt;
                opts.use_oldhash = input_parser.Results.use_oldhash;
                opts.use_oldcoeffmat = input_parser.Results.use_oldcoeffmat;
                opts.time_independent_A = input_parser.Results.time_independent_A;
                opts.time_independent_build_A = input_parser.Results.time_independent_build_A;
                opts.many_solution_capable = input_parser.Results.many_solution_capable;
                opts.max_lev = input_parser.Results.max_lev;
                opts.max_lev_coeffs = input_parser.Results.max_lev_coeffs;
                opts.adapt_threshold = input_parser.Results.adapt_threshold;
                opts.refinement_method = input_parser.Results.refinement_method;
                opts.adapt_initial_condition = input_parser.Results.adapt_initial_condition;
                opts.output_grid = input_parser.Results.output_grid;
                opts.save_output = input_parser.Results.save_output;
                opts.output_filename_id = input_parser.Results.output_filename_id;
                opts.plot_freq = input_parser.Results.plot_freq;
                opts.save_freq = input_parser.Results.save_freq;
                opts.use_connectivity = input_parser.Results.use_connectivity;
                opts.use_sparse_A = input_parser.Results.use_sparse_A;
                opts.case_ = input_parser.Results.case;
                opts.calculate_mass = input_parser.Results.calculate_mass;
                opts.normalize_by_mass = input_parser.Results.normalize_by_mass;
                opts.num_steps = input_parser.Results.num_steps;
                opts.start_time = input_parser.Results.start_time;
                
                opts.build_A = false;
                if sum(strcmp(opts.timestep_method,{'BE','CN','time_independent'}))>0
                    opts.build_A = true;
                end
                
                if opts.adapt
                    if opts.time_independent_build_A
                        disp("WARNING: setting 'time_independent_build_A',false to be compatible with 'adapt',true");
                    end
                    opts.time_independent_build_A = false;
                    if opts.time_independent_A
                        disp("WARNING: setting 'time_independent_A',false to be compatible with 'adapt',true");
                    end
                    opts.time_independent_A = false;
                end
                
                if opts.use_connectivity
                    if ~opts.use_oldhash
                        error("ERROR - must set 'use_oldhash' to use connectivity");
                    end
                    if opts.adapt
                        error("ERROR - cannot use adaptivity with use_connectivity==true");
                    end
                end
                
                if opts.normalize_by_mass
                    opts.calculate_mass = true;
                end
                
                if opts.adapt
                    opts.use_oldhash = false;
                    opts.calculate_mass = false;
                    opts.normalize_by_mass = false;
                end
                
                if ~strcmp(opts.output_grid,'quadrature')
                    if opts.normalize_by_mass
                        error("ERROR - cannot yet normalize by mass on output_grid other than 'quadrature'");
                    end
                    opts.calculate_mass = false;
                end
                
                if opts.max_lev_coeffs
                    if opts.use_oldcoeffmat
                        error("ERROR - max level coeffs not implemented for old coeff function");
                    end
                end
                
                if ~isempty(fieldnames(input_parser.Unmatched))
                    disp('Extra inputs:')
                    disp(input_parser.Unmatched)
                    error('Unrecognised input.')
                end
                
            end
            
        end
        
    end
    
end
