classdef EQUATION < handle
    
    properties
        unknown;
        terms;
        type = 'evolution';
        LHS_func;
        LHS_term;
    end
    
    methods
        
        function equation = EQUATION( unknown, terms, type, LHS_func )
            
            assert(strcmp(type,'evolution') || strcmp(type,'closure'),...
                   "equation type must be 'evolution' or 'closure' " )
            
            equation.unknown      = unknown;
            equation.terms        = terms;
            equation.type         = type;
            equation.unknown.type = equation.type;
            
            if( isempty(LHS_func) )
                
                num_dims = numel( unknown.dimensions );
                
                equation.LHS_func = cell( 1, num_dims + 1 );
                
                for i = 1 : num_dims
                    
                    equation.LHS_func{i} = @(x,p,t,d) 0*x+1;
                    
                end
                
                equation.LHS_func{num_dims+1} = @(t,p) 0*t+1;
                
            else
                
                equation.LHS_func = LHS_func;
                
            end
            
            sd_terms = cell( 1, num_dims ); % Does not include temporal component;
            
            for i = 1 : num_dims
                
                sd_terms{i} = SD_TERM( {MASS( equation.LHS_func{i} )} );
                
            end
            
            md_term = MD_TERM( num_dims, sd_terms );
            
            %For now assume mass matrix is time independent
            time_dep = false;
            equation.LHS_term = TERM( unknown, {unknown}, {md_term}, true, time_dep );
            
        end
        
        function [] = update_terms( obj, opts, t )
            
            obj.update_LHS_terms( opts, t )
            
            obj.update_RHS_terms( opts, t )
            
        end
        
        function [] = update_LHS_terms( obj, opts, t )
            
            obj.LHS_term.evaluate_coefficient_matrices( opts, t )
            
        end
        
        function [] = update_RHS_terms( obj, opts, t )
            
            for i = 1 : numel( obj.terms )
                
                term = obj.terms{i};
                
                if( term.linear )
                    
                    term.evaluate_coefficient_matrices( opts, t )
                    
                end
                
            end
            
        end
        
        function [ invMfval ] = MultiplyInverseMassMatrix( obj, opts, fval, t )
            
            Mx = @(x) obj.LHS_term.driver( opts, {x}, t );
            
            [ invMfval, ~, relres ] = pcg( Mx, fval, 1e-10, numel( fval ), [], [], fval );
            
            assert( relres < 1e-9, 'MultiplyInverseMassMatrix: pcg failed' )
            
        end
        
    end
end

