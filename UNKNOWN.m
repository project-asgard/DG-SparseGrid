classdef UNKNOWN < handle
    
    properties
        dimensions = {};
        deg;
        type = 'evolution';
        params = {};
        analytic_solutions = {};
        initial_conditions = {};
        hash_table;
        FMWT;
        transform_blocks;
        lo_global; % Start index of unknown in global solution vector
        hi_global; % End   index of unknown in global solution vector
    end
    
    methods
        
        function [ sz ] = size( obj )
            
            num_dims = numel(obj.dimensions);
            
            num_elements = numel(obj.hash_table.elements_idx);
            
            sz = obj.deg^num_dims * num_elements;
            
        end
        
        function fval_initial = get_initial_conditions( obj, opts, t )
            
            if ~exist( 't', 'var' )
                t = 0.0;
            end
            
            fval_initial...
                = md_eval_function( opts, opts.deg, obj.dimensions,[],...
                                    obj.initial_conditions, obj.hash_table,...
                                    obj.transform_blocks, t );
            
        end
        
        function fval_error = compute_error( obj, opts, fval_N, t )
            
            fval_A = md_eval_function( opts, opts.deg, obj.dimensions,[],...
                                       obj.analytic_solutions, obj.hash_table,...
                                       obj.transform_blocks, t );
            
            fval_error = sqrt( mean( ( fval_A(:) - fval_N(:) ).^2 ) );
            
        end
        
        function unknown = UNKNOWN( opts, dimensions, analytic_solutions, initial_conditions )
            
            unknown.dimensions         = dimensions;
            unknown.deg                = opts.deg;
            unknown.analytic_solutions = analytic_solutions;
            unknown.initial_conditions = initial_conditions;
            
            for i=1:numel(unknown.dimensions)
                unknown.dimensions{i}.lev = opts.lev(i);
            end
            
            [ unknown.hash_table.elements, unknown.hash_table.elements_idx ]...
                = hash_table_sparse_nD( unknown.get_lev_vec, opts.max_lev, opts.grid_type );
            
            [ ~, unknown.transform_blocks ] = OperatorTwoScale_wavelet2( opts.deg, opts.max_lev ); % Why max_lev here?
            
            switch numel(unknown.dimensions)
                case 1
                    [ unknown.FMWT, ~ ] = OperatorTwoScale_wavelet2( unknown.deg, unknown.dimensions{1}.lev );
                case 2
                    [ FMWT_x, ~ ]       = OperatorTwoScale_wavelet2( unknown.deg, unknown.dimensions{1}.lev );
                    [ FMWT_v, ~ ]       = OperatorTwoScale_wavelet2( unknown.deg, unknown.dimensions{2}.lev );
                    unknown.FMWT        = kron( FMWT_x, FMWT_v );
                otherwise
                    assert(false,'UKNOWN implemented for dimension<3')
            end
            
            compute_dimension_mass_mat( opts, unknown );
            
        end
        
        function lev_vec = get_lev_vec( obj )
            num_dims = numel(obj.dimensions);
            lev_vec  = zeros(num_dims,1);
            for d=1:num_dims
                lev_vec(d,1) = obj.dimensions{d}.lev;
            end
        end
        
        function set_bounds( obj, lo, hi )
            obj.lo_global = lo;
            obj.hi_global = hi;
        end
        
        function f_rs = convert_to_realspace( obj, f_w )
            
            f_rs = obj.FMWT' * f_w;
            
        end
        
        function f_w = convert_to_wavelet( obj, f_rs )
            
            f_w = obj.FMWT * f_rs;
            
        end
        
    end
end

