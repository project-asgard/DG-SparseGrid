classdef DIV < DIV_OR_GRAD
    
    methods
        
        function obj = DIV(varargin)   
            obj@DIV_OR_GRAD('div',varargin{:})   
            
            IBCL_ = obj.BCL;
            IBCR_ = obj.BCR;
            
            obj.IBCL = IBCL_;
            obj.IBCR = IBCR_;
        end
        
    end
    
end
