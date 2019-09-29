function [new_elem_lev_vecs, new_elem_pos_vecs, cnt] = ...
    get_element_children (lev_vec, pos_vec, max_lev, refinement_method)

%%
% Takes lev and pos vector as input, returns the same for a list of
% daughter elements which also obey one of the sparse adaption rules.

num_dimensions = numel(lev_vec);

assert (num_dimensions == numel (pos_vec) );

debug = 0;

new_elem_lev_vecs = [];
new_elem_pos_vecs = [];

%%
% Available refinement methods
%
% In this 2D example, assume we want to refine the (2,b) element and do so
% in all directions (dimensions)
%
%         1                a
%        / \              / \
%       2   3            b   c
%      / \              / \
%     4   5            d   e
%
% Method 1 (david)
% 
% number of new elements = 2 * num_dimension
%
% (4,b)
% (5,b)
% (2,d)
% (2,e)
%
% Method 2 (david + lin)
%
% number of new elements = 2 * num_dimension + 2^num_dimensions
%
% (4,b)
% (5,b)
% (2,d)
% (2,e)
% (4,d)
% (4,e)
% (5,d)
% (5,e)

if ~exist('method','var') || isempty(refinement_method)
    refinement_method = 2;
end

cnt = 0;
if refinement_method == 1 || refinement_method == 2
    
    for d=1:num_dimensions                   
                    
        %%
        % First daughter
        
        new_elem_lev_vec = lev_vec;
        new_elem_pos_vec = pos_vec;
                
        if new_elem_lev_vec(d)+1 <= max_lev
            
            new_elem_lev_vec(d) = new_elem_lev_vec(d)+1;
            new_elem_pos_vec(d) = new_elem_pos_vec(d)*2; % Assumes pos starts at 0
            
            new_elem_lev_vecs(cnt+1,:) = new_elem_lev_vec;
            new_elem_pos_vecs(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
            
            assert(new_elem_pos_vecs(cnt+1,d) >= 0);
            assert(new_elem_lev_vecs(cnt+1,d) >= 0);
            
            cnt = cnt + 1;
            
            new_lev_1D{d} = new_elem_lev_vec(d);
            new_pos_1D{d} = new_elem_pos_vec(d);
            
        end
        
        
        %%
        % If level is 1 or deeper add a second daughter
        
        if lev_vec(d) >= 1
                        
            new_elem_lev_vec = lev_vec;
            new_elem_pos_vec = pos_vec;
            
            if new_elem_lev_vec(d)+1 <= max_lev
                
                new_elem_lev_vec(d) = new_elem_lev_vec(d)+1;
                new_elem_pos_vec(d) = new_elem_pos_vec(d)*2+1; % Assumes pos starts at 0
                
                new_elem_lev_vecs(cnt+1,:) = new_elem_lev_vec;
                new_elem_pos_vecs(cnt+1,:) = new_elem_pos_vec; % Assumes pos starts at 0
                
                assert(new_elem_pos_vecs(cnt+1,d) >= 0);
                assert(new_elem_lev_vecs(cnt+1,d) >= 0);
                
                cnt = cnt + 1;
                
                new_lev_1D{d} = [new_lev_1D{d},new_elem_lev_vec(d)];
                new_pos_1D{d} = [new_pos_1D{d},new_elem_pos_vec(d)];
                
            end
            
        end
        
    end
    
end

if refinement_method == 2
    
    if num_dimensions == 1
        
        %%
        % Do nothing as method 2 == method 1 for 1D
        
    end
    
    % TODO : add max_lev check to method 2
    if num_dimensions == 2
        
        for i=1:numel(new_lev_1D{1})
            for j=1:numel(new_lev_1D{2})
                
                tmp_lev(1) = new_lev_1D{1}(i);
                tmp_lev(2) = new_lev_1D{2}(j);
                new_elem_lev_vecs(cnt+1,1:2) = tmp_lev;
                
                tmp_pos(1) = new_pos_1D{1}(i);
                tmp_pos(2) = new_pos_1D{2}(j);
                new_elem_pos_vecs(cnt+1,1:2) = tmp_pos;
                
                cnt = cnt + 1;
                
            end
        end
        
    end
    
    if num_dimensions == 3
        
        for i=1:numel(new_lev_1D{1})
            for j=1:numel(new_lev_1D{2})
                for k=1:numel(new_lev_1D{3})
                    
                    
                    tmp_lev(1) = new_lev_1D{1}(i);
                    tmp_lev(2) = new_lev_1D{2}(j);
                    tmp_lev(3) = new_lev_1D{3}(k);
                    
                    new_elem_lev_vecs(cnt+1,1:3) = tmp_lev;
                    
                    tmp_pos(1) = new_pos_1D{1}(i);
                    tmp_pos(2) = new_pos_1D{2}(j);
                    tmp_pos(3) = new_pos_1D{3}(k);
                    
                    new_elem_pos_vecs(cnt+1,1:3) = tmp_pos;
                    
                    cnt = cnt + 1;
                    
                end
            end
        end
        
    end
    
    if num_dimensions > 3
        
        %%
        % TODO
        
        assert(1==2);
        
    end
    
end

if cnt>0
    assert(min(new_elem_lev_vecs(:))>=0);
    assert(min(new_elem_pos_vecs(:))>=0);
end

end