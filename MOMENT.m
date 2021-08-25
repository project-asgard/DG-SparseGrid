classdef MOMENT
    
    %Source function used to capture moment integrals (u,g)_\W 
    %For now we're assuming the moment function g does not depend on time
    
    properties
        md_funcs
        fList
        vector
        moment_fval_integral
        moment_analytic_integral
    end
      
    methods
        function t = MOMENT(md_funcs_)
            t.md_funcs = md_funcs_;
            t.fList = {};
            t.vector = [];
            t.moment_fval_integral = [];
            t.moment_analytic_integral = [];
        end
        
        function [obj] = createFlist(obj,pde,opts)
            %Creates the coefficients of the moment vector on each domain.
            %No mass matrix inversion is done. 
            
            num_md_funcs = numel(obj.md_funcs);
            
            dims = pde.dimensions;
            num_dims  = numel(dims);
            
            obj.fList = cell(num_md_funcs,num_dims);
            
            for s=1:num_md_funcs
                md_func = obj.md_funcs{s};
                 
                for d=1:num_dims
                                       
                    construction_level = dims{d}.lev;
                    if opts.max_lev_coeffs
                        construction_level = opts.max_lev;
                    end
                    obj.fList{s,d} = forward_wavelet_transform(opts.deg,...
                        construction_level,dims{d}.min,dims{d}.max, ...
                        md_func{d},dims{d}.moment_dV,pde.params, ...
                        pde.transform_blocks, 0);                    
                    
                end
            end
        end
        
        function [obj] = createMomentVector(obj,opts,hash_table)
            %Actually contstructs the moment vector using fList.
            
            %Calculate only if adapt is true or the vector field is empty
            if isempty(obj.vector) || opts.adapt
                obj.vector = 0;
                num_md_funcs = numel(obj.md_funcs);

                for s=1:num_md_funcs
                    obj.vector = obj.vector + ...
                        combine_dimensions_D(opts.deg, obj.fList(s,:), 1, hash_table, opts.use_oldhash);
                end           
            end
            
        end
    end
    
end
