classdef GRAD < DIV_OR_GRAD
    
    methods
        
        function obj = GRAD(varargin)   
            obj@DIV_OR_GRAD('grad',varargin{:})     
            
            %Since we want the grad matrix to be a negative transpose of a
            %DIV matrix, we need to swap the wind direction as well as swap
            %the BCs N<=>D.  However, this swap will affect the BC call.
            %Instead we have another BC flag IBCL/IBCR which will build the
            %bilinear form with respect to Dirichlet/Free boundary
            %conditions will leaving the BC routine unaffected. 
            IBCL_ = obj.BCL;
            IBCR_ = obj.BCR;
            if strcmp(obj.BCL,'D')
                IBCL_ = 'N';
            end
            if strcmp(obj.BCL,'N')
                IBCL_ = 'D';
            end
            
            if strcmp(obj.BCR,'D')
                IBCR_ = 'N';
            end
            if strcmp(obj.BCR,'N')
                IBCR_ = 'D';
            end
            
            %Switch upwinding direction
            obj.LF = -obj.LF;
            
            obj.IBCL = IBCL_;
            obj.IBCR = IBCR_;
            
            
        end
        
    end
    
end
