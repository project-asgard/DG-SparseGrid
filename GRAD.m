classdef GRAD < PARTIAL_TERM_1D
    
    properties
        LF
        BCL
        BCR
        BCL_fList
        BCR_fList
    end
    
    methods
        function obj = GRAD(num_dims,g_,LF_,BCL_,BCR_,BCL_fList_,BCR_flist_,dat_)
            assert(nargin>=1);
            if nargin<2
                g_ = @(x,p,t,dat) x.*0+1;
            end
            if nargin<3
                LF_ = 0;
            end
            if nargin<4
                BCL_ = 'N';
            end
            if nargin<5
                BCR_ = 'N';
            end
            if nargin<6
                for d=1:num_dims % BC variation in all dimensions
                    BCL_fList_{d} = @(x,p,t) x.*0;
                end
                BCL_fList_{num_dims+1} = @(t,p) 1;  % time variation  
            end
            if nargin<7
                for d=1:num_dims % BC variation in all dimensions
                    BCR_flist_{d} = @(x,p,t) x.*0;
                end
                BCR_flist_{num_dims+1} = @(t,p) 1; 
            end
            if nargin<8
                dat_ = [];            
            end
            
            obj@PARTIAL_TERM_1D('grad',g_,dat_)
            
            obj.LF = LF_;
            obj.BCL = BCL_;
            obj.BCR = BCR_;
            obj.BCL_fList = BCL_fList_;
            obj.BCR_fList = BCR_flist_;
        end
    end
end