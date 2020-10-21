classdef TERM_ND
    
    properties
        name
        terms_1D
    end
    
    methods
        function t = TERM_ND(num_dims,terms_1D_)
            assert(nargin>0);
            if nargin<2
                for d=1:num_dims
                    terms_1D_{d} = TERM_1D();
                end
            end
            
            assert(num_dims==numel(terms_1D_));
            
            % Fill in empty dimensions with identity mass           
            for d=1:num_dims               
                if isempty(terms_1D_{d})                  
                    terms_1D_{d} = TERM_1D({MASS()});                 
                end             
            end
            
            t.terms_1D = terms_1D_;
            
        end
    end
    
end
