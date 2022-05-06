classdef TERM < handle
    
    properties
        linear = true;
        IMEX = 'EX';
        time_dependent = true;
        output_unknown;
        input_unknowns;
        input_g2l; % --- global to local map for unknowns input to term
        descriptor;
        A_data;
    end
    
    methods
        
        function [ F ] = driver( obj, opts, t, x )
            
            if( nargin < 4 )
                apply_to_x = false;
            else
                apply_to_x = true;
                for i = 1 : numel( obj.input_unknowns )
                    assert( all(size(x{i}) == size(obj.input_unknowns{i}.fval)) )
                end
            end
            
            % --- Called by Time-Stepper ---
            
            F = zeros(size(obj.output_unknown.fval));
            
            for i = 1 : numel(obj.input_unknowns)
                
                obj.A_data{i} = matrix_assembly_data( obj.input_unknowns{i}, obj.output_unknown, opts );
                
                if( isa( obj.descriptor{i}, 'MD_TERM' ) )
                        
                    if( obj.time_dependent )

                        obj.descriptor{i}...
                            = obj.get_coeff_MATS( opts, obj.descriptor{i}, t ); % Needs to have same cell structure as input unknowns

                        % Generates 1D stiffness matrix for each dimension
                        %   Assumes input and output variables have the
                        %   same dimensionality.
                        if( apply_to_x )
                            F = F + apply_A_term( opts, obj.descriptor{i}, obj.A_data{i}, x{i}                      , obj.output_unknown.deg );
                        else
                            F = F + apply_A_term( opts, obj.descriptor{i}, obj.A_data{i}, obj.input_unknowns{i}.fval, obj.output_unknown.deg );
                        end

                    else % --- time independent

                    end
                        
                end
            end
            
        end
        
        function [ md_term ] = get_coeff_MATS( obj, opts, md_term, t ) % Ignore any LHS stuff (handeled in EQUATION class)
            
            num_dims = numel( obj.output_unknown.dimensions );
            for d = 1:num_dims
                
                dim = obj.output_unknown.dimensions{d};
                term_1D = md_term.terms_1D{d};

                construction_level = opts.lev;
                if opts.max_lev_coeffs && ~term_1D.time_dependent
                    construction_level = opts.max_lev;
                end
                
                md_term.terms_1D{d}...
                    = obj.sd_term_coeff_matrix( opts.deg, t, dim, term_1D, obj.output_unknown.params, ...
                        obj.output_unknown.transform_blocks, construction_level );
                
            end
            
        end
        
        function [ term_1D ] = sd_term_coeff_matrix( obj, deg, t, dim, term_1D, params, FMWT_blocks, construction_level )
            
            % Get the matrix for each pterm

            for i=1:numel(term_1D.pterms)
                term_1D.pterms{i} = obj.pterm_coeff_matrix(deg,t,dim,term_1D.pterms{i},params,FMWT_blocks,construction_level,logical(i-1));
            end

            % Chain together the partial terms

            mat = eye(size(term_1D.pterms{1}.mat));
            
            for i=1:numel(term_1D.pterms)
                mat = mat * term_1D.pterms{i}.mat;
            end

            term_1D.mat = mat;
            
        end
        
        function [pterm] = pterm_coeff_matrix(obj,deg,t,dim,pterm,params, FMWT_blocks, coeff_level, apply_mass_inverse )

            % Get the RHS matrix (L) for this pterm

            [L,~] = coeff_matrix(deg,t,dim,pterm,params,FMWT_blocks,coeff_level); % note that dV is applied in here

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

                pterm.mat = pterm.LHS_mass_mat\L;
                
            else
                
                pterm.mat = L;
            
            end

            % Store M for use in application to the BC

        end
        
        function term = TERM( output_unknown, input_unknowns, descriptor, linear, time_dependent, IMEX )
            
            term.output_unknown = output_unknown;
            term.input_unknowns = input_unknowns;
            term.descriptor     = descriptor;
            
            if exist( 'linear', 'var' )
                assert( (linear==true | linear==false), 'linear must be boolean' );
                term.linear = linear;
            end
            
            if exist( 'time_dependent', 'var' )
                assert( (time_dependent==true | time_dependent==false), 'time_dependent must be boolean' );
                term.time_dependent = time_dependent;
            end
            
            if exist( 'IMEX', 'var' )
                term.IMEX = IMEX;
            end
            
        end
        
    end
end

