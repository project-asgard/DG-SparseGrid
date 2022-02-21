classdef MD_TERM
    
    properties
        name
        terms_1D
        imex
    end
    
    methods
        function t = MD_TERM(num_dims,terms_1D_,imex_)
            assert(nargin>0);
            if nargin<2
                for d=1:num_dims
                    terms_1D_{d} = SD_TERM();
                end
                imex_ = 'N';
            end
            
            assert(num_dims==numel(terms_1D_));
            
            % Fill in empty dimensions with identity mass           
            for d=1:num_dims               
                if isempty(terms_1D_{d})                  
                    terms_1D_{d} = SD_TERM({MASS()});                 
                end             
            end
            
            t.terms_1D = terms_1D_;
            t.imex = imex_;
            
            
        end
    end
    
end
