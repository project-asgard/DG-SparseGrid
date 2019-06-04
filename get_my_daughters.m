function [new_elem_lev_vecs, new_elem_pos_vecs, cnt] = get_my_daughters (lev_vec, pos_vec, max_lev)

%%
% Takes lev and pos vector as input, returns the same for a list of
% daughter elements which also obey the sparse selecttion rule for the next
% deepest level.

num_dimensions = numel(lev_vec);

assert (num_dimensions == numel (pos_vec) );

debug = 0;

new_elem_lev_vecs = [];
new_elem_pos_vecs = [];

%%
% Methods
%
% 1 = add only 2 elements per dimension

method = 1;

if method == 1
    
    cnt = 0;
    for d=1:num_dimensions
                    
                    
        %%
        % First daughter
        
        new_elem_lev_vec = lev_vec;
        new_elem_pos_vec = pos_vec;
        new_elem_lev_vec(d) = new_elem_lev_vec(d)+1;
        new_elem_pos_vec(d) = new_elem_pos_vec(d)*2; % Assumes pos starts at 0
        
        new_elem_lev_vecs(cnt+1,:) = new_elem_lev_vec;
        new_elem_pos_vecs(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
        
        assert(new_elem_pos_vecs(cnt+1,d) >= 0);
        assert(new_elem_lev_vecs(cnt+1,d) >= 0);
        
        cnt = cnt + 1;
        
        %%
        % If level is 1 or deeper add a second daughter
        
        if lev_vec(d) >= 1
            
            new_elem_lev_vec = lev_vec;
            new_elem_pos_vec = pos_vec;
            
            new_elem_lev_vec(d) = new_elem_lev_vec(d)+1;
            new_elem_pos_vec(d) = new_elem_pos_vec(d)*2+1; % Assumes pos starts at 0
               
            new_elem_lev_vecs(cnt+1,:) = new_elem_lev_vec;
            new_elem_pos_vecs(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
            
            assert(new_elem_pos_vecs(cnt+1,d) >= 0);
            assert(new_elem_lev_vecs(cnt+1,d) >= 0);
            
            cnt = cnt + 1;
            
        end
        
    end
    
end

end