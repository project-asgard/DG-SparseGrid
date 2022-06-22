function [pde,hash_table,fval,A_data] ...
    = adapt_stripped(pde,opts,hash_table,fval,coarsen_or_refine)

num_elements    = numel(hash_table.elements_idx);
num_dims  = numel(pde.dimensions);
deg             = opts.deg;

% coarsen_ = 0;
% refine_ = 0;

element_DOF = deg^num_dims;

relative_threshold = opts.adapt_threshold;

refine_threshold  = max(abs(fval)) * relative_threshold;
coarsen_threshold = refine_threshold * 0.1;

new_element_value = 1e-15;

%Get PDE lev vec so we dont go over it
for d=1:num_dims
    pde_lev_vec(d) = pde.dimensions{d}.lev;
end

debug   = 0;
refinement_method  = opts.refinement_method; % 1 = david, 2 = lin; see get_child_elements for description

coarsen = false;
refine  = false;
switch coarsen_or_refine
    case 'r'
        refine = true;
    case 'c'
        coarsen = true;
    otherwise
        error("adapt algorithm must be told to refine 'r' or coarsen 'c'.  Please specify.")
end
        

if refine_threshold <= 0 % just don't do anything for zero valued functions
    coarsen = 0;
    refine = 0;
    refine_previous = 0;
end

assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));

%%
% Coarsen

if coarsen
    
    if ~opts.quiet
        disp('Coarsening ...')
        fprintf('    Initial number of elements: %i\n', numel(hash_table.elements_idx));
        fprintf('    Initial number of DOFs: %i\n', numel(hash_table.elements_idx)*deg^num_dims);
    end
    
    num_remove = 0;
    elements_to_remove = [];
    
    num_new_leaf_elements = 0;
    new_leaf_elements = [];
    
    for n=1:num_elements
        
        idx = hash_table.elements_idx(n);
        lev_vec = hash_table.elements.lev_p1(idx,:)-1;
        pos_vec = hash_table.elements.pos_p1(idx,:)-1;
        
        assert(max(lev_vec)<=opts.max_lev);
        
        [lev_vec_, pos_vec_] = md_idx_to_lev_pos (num_dims, opts.max_lev, double(idx)); %IDK why the double() is required.  Hack for now!!!
        assert(norm(lev_vec-lev_vec_)==0);
        assert(norm(pos_vec-pos_vec_)==0);
        
        gidx1 = (n-1)*element_DOF+1;
        gidx2 = n*element_DOF;
        
        element_sum = sqrt(sum(fval(gidx1:gidx2).^2));
        element_max = max(abs(fval(gidx1:gidx2)));
        
        %%
        % check if the element needs refining, if it is at least level 1,
        % and is labeled as a leaf
        
        if element_max <= coarsen_threshold ...
                && min(hash_table.elements.lev_p1(idx,:)>=1) % level must be >= 0 at present
            %&& hash_table.elements.type(idx) == 2
            
            %%
            % get element children and check if any are live elements
            
            %[num_live_children, has_complete_children] = ...
            %    number_of_live_children (hash_table, lev_vec, pos_vec, opts.max_lev, refinement_method);
            
            %%
            % only coarsen (remove) this element if it has no (live)
            % daughters
            
            %if num_live_children == 0
            
            if debug
                disp(['    Removing : ',num2str(hash_table.elements.lev_p1(idx,:)-1)]);
                disp(['        its type is : ', num2str(hash_table.elements.type(idx))]);
            end
            
            num_remove = num_remove + 1;
            elements_to_remove(num_remove) = n;
            
            %%
            % determine level above leaf nodes and label them
            
            %                 parent_elements_idx = get_parent_elements_idx(hash_table, idx, opts.max_lev, refinement_method );
            %
            %                 for ii=1:numel(parent_elements_idx)
            %
            %                     % make sure the element we want to be a leaf is already in
            %                     % the table and active
            %                     assert(hash_table.elements.type( parent_elements_idx(ii) ) >= 1);
            %
            %                     % store the elements which will become leafs below
            %                     num_new_leaf_elements = num_new_leaf_elements + 1;
            %                     new_leaf_elements(num_new_leaf_elements) = parent_elements_idx(ii);
            %
            %                 end
            
            %end
            
        end
    end
    
    
    
    %%
    % Now remove elements
    
    assert(numel(elements_to_remove)==num_remove);
    
    remove_DOF_list = [];
    for n=1:num_remove
        
        %%
        % Remove entries from element table (recall sparse storage means =0
        % removes it from the table
        
        hash_table.elements.lev_p1(hash_table.elements_idx(elements_to_remove(n)),:) = 0;
        hash_table.elements.pos_p1(hash_table.elements_idx(elements_to_remove(n)),:) = 0;
        hash_table.elements.type(hash_table.elements_idx(elements_to_remove(n)))= 0;
        
        %%
        % Add this elements DOF to the list to be removed from fval
        
        nn = elements_to_remove(n);
        
        i1 = (nn-1)*element_DOF+1; % Get the start and end global row indices of the element
        i2 = (nn)*element_DOF;
        assert(i2-i1==element_DOF-1);
        
        remove_DOF_list = [remove_DOF_list,i1:i2];
        
    end
    
    %%
    % Remove elements from elements_idx, and DOF from fval
    
    hash_table.elements_idx(elements_to_remove) = [];
    fval(remove_DOF_list) = [];
    
    if ~opts.quiet
        fprintf('    Final number of elements: %i\n', numel(hash_table.elements_idx));
        fprintf('    Final number of DOFs: %i\n', numel(hash_table.elements_idx)*deg^num_dims);
    end
    
    %%
    % Label new leaf elements
    
    %     for n=1:num_new_leaf_elements
    %         idx = new_leaf_elements(n);
    %
    %         %%
    %         % assert that the element we are making a leaf does not have a
    %         % complete set of live children
    %         lev_vec = hash_table.elements.lev_p1(idx,:)-1;
    %         pos_vec = hash_table.elements.pos_p1(idx,:)-1;
    %         [num_live_children, has_complete_children] = ...
    %             number_of_live_children (hash_table, lev_vec, pos_vec, opts.max_lev, refinement_method);
    %         if ~has_complete_children
    %             hash_table.elements.type(idx) = 2;
    %         end
    %     end
    
end

%%
% Refine

if refine
    
    if ~opts.quiet
        disp('Refining ...')
        fprintf('    Initial number of elements: %i\n', numel(hash_table.elements_idx));
        fprintf('    Initial number of DOFs: %i\n', numel(hash_table.elements_idx)*deg^num_dims);
    end
    
    num_elements = numel(hash_table.elements_idx);
    cnt = 0;
    clear new_elements_lev_vec;
    clear new_elements_pos_vec;
    for n=1:num_elements
        
        idx = hash_table.elements_idx(n);
        
        gidx1 = (n-1)*element_DOF+1;
        gidx2 = n*element_DOF;
        
        element_sum = sqrt(sum(fval(gidx1:gidx2).^2));
        element_max = max(abs(fval(gidx1:gidx2)));
        
        %%
        % Check for refinement
        
        if element_max >= refine_threshold %&& hash_table.elements.type(idx) == 2
            
            if debug; disp([...
                    '    refine ? yes, fval = ', num2str(element_max,'%1.1e'), ...
                    ', type = ', num2str(hash_table.elements.type(idx)), ...
                    ', lev_vec = ', num2str(hash_table.elements.lev_p1(idx,:)-1) ...
                    ', pos_vec = ', num2str(hash_table.elements.pos_p1(idx,:)-1) ...
                    ', idx = ', num2str(idx) ...
                    ]); end
            
            [child_elements_idx, num_children] = ...
                get_child_elements_idx(num_dims, opts.max_lev, idx, refinement_method);
            
            if num_children > 0
                
                if debug
                    for nn=1:num_children
                        [lev_vec, pos_vec] = md_idx_to_lev_pos(num_dims, opts.max_lev, child_elements_idx(nn));
                        disp(['        adding element with lev : ',num2str(lev_vec), ...
                            ', idx = ', num2str(child_elements_idx(nn))]);
                    end
                end
                
                new_elements_idx(cnt+1:cnt+num_children) = child_elements_idx;
                hash_table.elements.type(idx) = 1; % Now that this element has been refined it is no longer a leaf.
                cnt = cnt + num_children;
                
            end
            
        else
            
            if debug; disp(['    refine ?  no, fval = ', num2str(element_max,'%1.1e'), ...
                    ' type = ', num2str(hash_table.elements.type(idx))]); end
            
        end
        
    end
    
    %%
    % Now add these elements with (almost) zero coefficient to the
    % elements table and elementsIDX
    
    num_try_to_add = cnt;
    num_elements_added = 0;
    for i=1:num_try_to_add
        
        idx = new_elements_idx(i);
        
        %%
        % Sanity check the new element
        
        assert(idx>=0);
        
        %%
        % If element does not exist, add its idx to the list of active elements
        % (hash_table.elements_idx)
        
        [lev_vec, pos_vec] = md_idx_to_lev_pos(num_dims, opts.max_lev, idx);
        if hash_table.elements.type(idx) == 0 && all(lev_vec <= pde_lev_vec) % element not already enabled and level does not grow
            
            num_elements_added = num_elements_added + 1;
            position_in_elements_idx = num_elements+num_elements_added;
            hash_table.elements_idx(position_in_elements_idx) = idx; % Extend element list
            i1 = (position_in_elements_idx-1)*element_DOF+1; % Get the start and end global row indices of the new element
            i2 = (position_in_elements_idx)*element_DOF;
            assert(i2-i1==element_DOF-1);
            fval(i1:i2) = new_element_value; % Extend coefficient list with near zero magnitude (ideally would be zero)
            
            %[lev_vec, pos_vec] = md_idx_to_lev_pos(num_dims, opts.max_lev, idx);
            
            hash_table.elements.lev_p1(idx,:) = lev_vec+1; % NOTE : have to start lev  index from 1 for sparse storage
            hash_table.elements.pos_p1(idx,:) = pos_vec+1; % NOTE : have to start cell index from 1 for sparse storage
            hash_table.elements.type(idx) = 1;
            
        end
        
        %         %%
        %         % Set element to type to leaf
        %
        %         [num_live_children, has_complete_children] = ...
        %             number_of_live_children_idx (hash_table, idx, opts.max_lev, refinement_method);
        %
        %         if has_complete_children
        %             hash_table.elements.type(idx) = 1;
        %         else
        %             hash_table.elements.type(idx) = 2;
        %         end
        
    end
    
    %     if ~opts.quiet; fprintf('    Refine on : added %i elements\n', num_elements_added); end
    
    if ~opts.quiet
        fprintf('    Final number of elements: %i\n', numel(hash_table.elements_idx));
        fprintf('    Final number of DOFs: %i\n', numel(hash_table.elements_idx)*deg^num_dims);
    end
    
end

%%
% Some more sanity checking on the refined element list / fval

assert(numel(fval)==numel(hash_table.elements_idx)*element_DOF);
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));

for i=1:numel(hash_table.elements_idx)
    lev_vec = hash_table.elements.lev_p1(hash_table.elements_idx(i));
    assert(min(lev_vec>0));
end


elements_idx0 = hash_table.elements_idx;

%%
% Update all the setup outputs which need updating on the new element list


%% Update dims and coeffs
% Update the time-indepedent coeff mats to the new size
lev_vec = zeros(num_dims, 1);
for d=1:num_dims
    lev_vec(d) = max(hash_table.elements.lev_p1(:,d)-1);
end

%pde = compute_dimension_mass_mat(opts,pde);

if opts.max_lev_coeffs
    pde = get_coeff_mats_rechain(pde, deg, lev_vec);
end

for d=1:num_dims
    pde.dimensions{d}.lev = lev_vec(d);
end

assert(norm(pde.get_lev_vec-lev_vec)==0);

% If we don't want to store the max lev coeffs, regen them
if ~opts.max_lev_coeffs
    t = 0;
    TD = 0;
    pde = compute_dimension_mass_mat(opts,pde);
    pde = get_coeff_mats(pde,opts,t,TD);
end

%% Update A_data

A_data = global_matrix(pde,opts,hash_table);


%% More checks
assert(numel(fval)==numel(hash_table.elements_idx)*element_DOF);
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));
assert(sum(hash_table.elements_idx-elements_idx0)==0);

end
