classdef DIV_OR_GRAD < PARTIAL_SD_TERM
    
    properties
        LF
        BCL
        BCR
        BCL_fList
        BCR_fList
        IBCL
        IBCR
    end
    
    methods
        function obj = DIV_OR_GRAD(type_,num_dims,g_,LHS_mass_g_,LF_,BCL_,BCR_,BCL_fList_,BCR_fList_,dat_,dV_)
            assert(nargin>=2);
            if nargin<3
                g_ = @(x,p,t,dat) x.*0+1;
            end
            
            LHS_mass_g = @(x,p,t,dat) x.*0+1;
            if exist('LHS_mass_g_','var')
                if ~isempty(LHS_mass_g_)
                    LHS_mass_g = LHS_mass_g_;
                end
            end
            
            
            if nargin<5
                LF_ = 0;
            end
            if nargin<6
                BCL_ = 'N';
            end
            if nargin<7
                BCR_ = 'N';
            end
            
            
            for d=1:num_dims % BC variation in all dimensions
                BCL_fList{d} = @(x,p,t) x.*0;
            end
            BCL_fList{num_dims+1} = @(t,p) 1;  % time variation
            if exist('BCL_fList_','var')
                if ~isempty(BCL_fList_)
                    BCL_fList = BCL_fList_;
                end
            end
            
            for d=1:num_dims % BC variation in all dimensions
                BCR_fList{d} = @(x,p,t) x.*0;
            end
            BCR_fList{num_dims+1} = @(t,p) 1;  % time variation
            if exist('BCR_fList_','var')
                if ~isempty(BCR_fList_)
                    BCR_fList = BCR_fList_;
                end
            end
            
            dat = [];
            if exist('dat_','var')
                if ~isempty(dat_)
                    dat = dat_;
                end
            end
                      
            
            if nargin<11
                dV_ = @(x,p,t,d) x.*0+1;
            end
            
            obj@PARTIAL_SD_TERM(type_,g_,LHS_mass_g,dat,dV_)
            
            obj.LF = LF_;
            obj.BCL = BCL_;
            obj.BCR = BCR_;
            obj.BCL_fList = BCL_fList;
            obj.BCR_fList = BCR_fList;
        end
    end
end
