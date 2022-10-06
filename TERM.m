classdef TERM < handle
    
    properties
        linear = true;
        IMEX = 'EX';
        time_dependent = true;
        output_unknown;
        input_unknowns;
        descriptor;
    end
    
    methods
        
        function obj = TERM( output_unknown, input_unknowns, descriptor, linear, time_dependent, IMEX )
            
            obj.output_unknown = output_unknown;
            obj.input_unknowns = input_unknowns;
            obj.descriptor     = descriptor;
            
            if exist( 'linear', 'var' ) && ~isempty( linear )
                assert( (linear==true | linear==false), 'linear must be boolean' );
                obj.linear = linear;
            end
            
            if exist( 'time_dependent', 'var' ) && ~isempty( time_dependent )
                assert( (time_dependent==true | time_dependent==false), 'time_dependent must be boolean' );
                obj.time_dependent = time_dependent;
            end
            
            if exist( 'IMEX', 'var' ) && ~isempty( IMEX )
                obj.IMEX = IMEX;
            end
            
        end
        
        function [ F ] = driver( obj, opts, Q, t )
            
            assert( numel(Q) == numel(obj.input_unknowns) )
            
            % --- Called by Time-Stepper ---
            
            F = zeros(obj.output_unknown.size(),1);
            
            if( obj.linear )
                
                for i = 1 : numel(obj.input_unknowns)
                    
                    if( isa( obj.descriptor{i}, 'MD_TERM' ) )
                    
                        F = F + apply_A_term( opts, obj.descriptor{i}, obj.input_unknowns{i}, obj.output_unknown, Q{i} );
                                              
                    end
                    
                end
                
            else
                
                Q_rs = cell(numel(Q),1);
                
                for i = 1 : numel(obj.input_unknowns)
                    
                    unknown = obj.input_unknowns{i};
                    
                    Q_rs{i} = unknown.convert_to_realspace( Q{i} );
                    
                end
                
                F_rs = obj.descriptor{1}( opts, Q_rs, t );
                
                F = obj.output_unknown.convert_to_wavelet( F_rs );
                
            end
            
        end
        
        function evaluate_coefficient_matrices( obj, opts, t )
            
            for i = 1 : numel(obj.input_unknowns)
                
                if( isa( obj.descriptor{i}, 'MD_TERM' ) )
                    
                    if( or( obj.time_dependent, isempty(obj.descriptor{i}.terms_1D{1}.mat) ) )
                        
                        obj.descriptor{i} = obj.get_coeff_MATS( opts, obj.descriptor{i}, obj.input_unknowns{i}, t ); % Needs to have same cell structure as input unknowns
                        
                    end
                
                end
                
            end
            
        end
        
        function [ md_term ] = get_coeff_MATS( obj, opts, md_term, input_unknown, t ) % Ignore any LHS stuff (handeled in EQUATION class)
            
            num_dims = numel( input_unknown.dimensions );
            for d = 1:num_dims
                
                dim = input_unknown.dimensions{d};
                term_1D = md_term.terms_1D{d};

                %construction_level = opts.lev; % --- Should be obj.output_unknown.dimensions{d}.lev?
                construction_level = input_unknown.dimensions{d}.lev;
                if opts.max_lev_coeffs && ~term_1D.time_dependent
                    construction_level = opts.max_lev;
                end
                
                md_term.terms_1D{d}...
                    = obj.sd_term_coeff_matrix( opts.deg, t, dim, term_1D, input_unknown.params, ...
                        input_unknown.transform_blocks, construction_level );
                
            end
            
        end
        
        function [ term_1D ] = sd_term_coeff_matrix( obj, deg, t, dim, term_1D, params, FMWT_blocks, construction_level )
            
            % Get the matrix for each pterm

            for i=1:numel(term_1D.pterms)
                term_1D.pterms{i} = obj.pterm_coeff_matrix(deg,t,dim,term_1D.pterms{i},params,FMWT_blocks,construction_level,logical(i-1));
            end

            % Chain together the partial terms

            mat = term_1D.pterms{1}.mat;
            
            for i=2:numel(term_1D.pterms)
                mat = mat * term_1D.pterms{i}.mat;
            end

            term_1D.mat = mat;
            
        end
        
        function [ pterm ] = pterm_coeff_matrix( obj, deg, t, dim, pterm, params, FMWT_blocks, coeff_level, apply_mass_inverse )

            % Get the RHS matrix (L) for this pterm

            if( strcmp( pterm.type, 'moment' ) )
                L = forward_wavelet_transform...% Store Transpose
                      ( deg, coeff_level, dim.min, dim.max, pterm.g, pterm.dV, params, FMWT_blocks, t )';
            else
                [ L, ~ ] = coeff_matrix( deg, t, dim, pterm, params, FMWT_blocks, coeff_level ); % note that dV is applied in here
            end

            % Get the LHS mass matrix (M) for this pterm
            % this mass matrix needs to have the volume jacobian incorporated instead
            % of the possibly different surface jacobian.
            if isempty(pterm.LHS_mass_mat)
                lhs_mass_pterm = MASS(pterm.LHS_mass_g,'','',dim.moment_dV);
                [M,~] = coeff_matrix(deg,t,dim,lhs_mass_pterm,params,FMWT_blocks,coeff_level); 
                pterm.LHS_mass_mat = full(M);
            end

            if( apply_mass_inverse )
                
                % Move M to the RHS

                pterm.mat = pterm.LHS_mass_mat \ L;
                
            else
                
                pterm.mat = L;
            
            end

            % Store M for use in application to the BC

        end
        
    end
end

