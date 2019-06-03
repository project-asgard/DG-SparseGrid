function [new_elem_lev_vecs, new_elem_pos_vecs, cnt] = getMyDaughters (lev_vec, pos_vec)

%%
% Takes lev and pos vector as input, returns the same for a list of
% daughter elements which also obey the sparse selecttion rule for the next
% deepest level. 

nDims = numel(lev_vec);

assert (nDims == numel (pos_vec) );

debug = 1;

if lev_vec(1) == 2
    disp('should never get here');
end

new_elem_lev_vecs = [];
new_elem_pos_vecs = [];

%%
% Generate list of new elements which satisfy selection rule

cnt = 0;
for d=1:nDims
    
    %%
    % Add two daughter nodes in each dimension for each leaf
    % element
    
    for d2=1:nDims
        
        %%
        % First daughter
          
        new_elem_lev_vec = lev_vec;
        new_elem_pos_vec = pos_vec;
        new_elem_lev_vec(d2) = new_elem_lev_vec(d)+1;
        new_elem_pos_vec(d2) = new_elem_pos_vec(d)*2; % Assumes pos starts at 0
        
        if sum(new_elem_lev_vec)<=lev_vec(d)+1 & new_elem_lev_vec(d2)<=pde.maxLev % Sparse grid selection rule AND max depth check
            
            new_elem_lev_vecs(cnt+1,:) = new_elem_lev_vec;
            new_elem_pos_vecs(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
            
            assert(new_elem_pos_vecs(cnt+1,d2) >= 0);
            assert(new_elem_lev_vecs(cnt+1,d2) >= 0);
            
            cnt = cnt + 1;
                        
        else
            
            disp('element not added because it did not obey sparse selection rule or is > pde.maxLev');
            
        end
        
        %%
        % Second daughter
        
        new_elem_lev_vec = lev_vec;
        new_elem_pos_vec = pos_vec;
        
        new_elem_lev_vec(d2) = new_elem_lev_vec(d)+1;
        new_elem_pos_vec(d2) = new_elem_pos_vec(d)*2+1; % Assumes pos starts at 0
        
        if sum(new_elem_lev_vec)<=lev_vec(d)+1 & new_elem_lev_vec(d2)<=pde.maxLev % Sparse grid selection rule AND max depth check
            
            new_elem_lev_vecs(cnt+1,:) = new_elem_lev_vec;
            new_elem_pos_vecs(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
            
            assert(new_elem_pos_vecs(cnt+1,d2) >= 0);
            assert(new_elem_lev_vecs(cnt+1,d2) >= 0);
            
            cnt = cnt + 1;
                        
        else
            
            disp('element not added because it did not obey sparse selection rule');
            
            
        end
        
    end
    
end

end