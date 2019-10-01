function [pde,fval,hash_table,A_data,Meval,nodes,coord,fval_previous] ...
    = adapt(pde,opts,fval,hash_table,nodes0,fval_realspace0,coarsen_,refine_,fval_previous)

num_elements    = numel(hash_table.elements_idx);
num_dims  = numel(pde.dimensions);
deg             = pde.deg;

element_DOF = deg^num_dims;

relative_threshold = 1e-3;

refine_threshold  = max(abs(fval)) * relative_threshold;
coarsen_threshold = refine_threshold * 0.1;

new_element_value = 1e-15;

debug   = 0;
refinement_method  = 1; % 1 = david, 2 = lin; see get_child_elements for description

coarsen = 0;
if exist('coarsen_','var') && ~isempty(coarsen_)
    coarsen = coarsen_;
end

refine  = 0;
if exist('refine_','var') && ~isempty(refine_)
    refine = refine_;
end

refine_previous = 0;
if exist('fval_previous','var') && ~isempty(fval_previous)
    refine_previous = 1;
    assert(numel(fval)==numel(fval_previous));
end

if refine_threshold <= 0 % just don't do anything for zero valued functions
    coarsen = 0;
    refine = 0;
    refine_previous = 0;
end

assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));

%%
% Store unrefined fval for comparison with refined fval after

pde0  = pde;
fval0 = fval;
for d=1:num_dims
    [Meval0{d},nodes0{d}] = matrix_plot_D(pde,pde.dimensions{d});
end
fval_realspace0 = wavelet_to_realspace(pde0,opts,Meval0,fval0,hash_table);

%%
% Plot the grid (1 and 2D only)

plot_grid = 1;
if plot_grid && ~opts.quiet
    plot_adapt(pde,opts,hash_table,1);
    plot_adapt_triangle(pde,opts,hash_table,7);
end

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
        
        gidx1 = (n-1)*element_DOF+1;
        gidx2 = n*element_DOF;
        
        element_sum = sqrt(sum(fval(gidx1:gidx2).^2));
        
        %%
        % check if the element needs refining, if it is at least level 1,
        % and is labeled as a leaf
        
        if element_sum <= coarsen_threshold ...
                && min(hash_table.elements.lev_p1(idx,:))>=2 ...
                && hash_table.elements.type(idx) == 2
            
            %%
            % get element children and check if any are live elements
            
            [num_live_children, has_complete_children] = ...
                number_of_live_children (hash_table, lev_vec, pos_vec, pde.max_lev, refinement_method);
            
            %%
            % only coarsen (remove) this element if it has no (live)
            % daughters
            
            if num_live_children == 0 
                
                if debug
                    disp(['    Removing : ',num2str(hash_table.elements.lev_p1(idx,:)-1)]);
                    disp(['        its type is : ', num2str(hash_table.elements.type(idx))]);
                end
                
                num_remove = num_remove + 1;
                elements_to_remove(num_remove) = n;
                
                %%
                % determine level above leaf nodes and label them
                
                parent_elements_idx = get_parent_elements_idx(hash_table, idx, pde.max_lev, refinement_method );
                
                for ii=1:numel(parent_elements_idx)
                    
                    % make sure the element we want to be a leaf is already in
                    % the table and active
                    assert(hash_table.elements.type( parent_elements_idx(ii) ) >= 1);
                    
                    % store the elements which will become leafs below
                    num_new_leaf_elements = num_new_leaf_elements + 1;                   
                    new_leaf_elements(num_new_leaf_elements) = parent_elements_idx(ii);
                    
                end
                
            end
            
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
    
    for n=1:num_new_leaf_elements
        idx = new_leaf_elements(n);
        
        %%
        % assert that the element we are making a leaf does not have a
        % complete set of live children
        lev_vec = hash_table.elements.lev_p1(idx,:)-1;
        pos_vec = hash_table.elements.pos_p1(idx,:)-1;
        [num_live_children, has_complete_children] = ...
            number_of_live_children (hash_table, lev_vec, pos_vec, pde.max_lev, refinement_method);
        if ~has_complete_children
            hash_table.elements.type(idx) = 2;
        end
    end
    
end

assert(numel(fval)==numel(hash_table.elements_idx)*element_DOF);
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));
num_elements = numel(hash_table.elements_idx);


%%
% Plot the coarsened grid (1 and 2D only)

plot_grid = 1;
if plot_grid && ~opts.quiet
    plot_adapt(pde,opts,hash_table,2);
    plot_adapt_triangle(pde,opts,hash_table,8);
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
        
        %%
        % Check for refinement
        
        if element_sum >= refine_threshold && hash_table.elements.type(idx) == 2
            
            if debug; disp([...
                    '    refine ? yes, fval = ', num2str(element_sum,'%1.1e'), ...
                    ', type = ', num2str(hash_table.elements.type(idx)), ...
                    ', lev_vec = ', num2str(hash_table.elements.lev_p1(idx,:)-1) ...
                    ', pos_vec = ', num2str(hash_table.elements.pos_p1(idx,:)-1) ...
                    ]); end
            
            [child_elements_idx, num_children] = ...
                get_child_elements_idx(hash_table, idx, pde.max_lev, refinement_method);
            
            if num_children > 0
                
                new_elements_idx(cnt+1:cnt+num_children) = child_elements_idx;              
                hash_table.elements.type(idx) = 1; % Now that this element has been refined it is no longer a leaf.               
                cnt = cnt + num_children;
                
            end
            
        else
            
            if debug; disp(['    refine ?  no, fval = ', num2str(element_sum,'%1.1e'), ...
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
        
        if hash_table.elements.type(idx) == 0
            
            num_elements_added = num_elements_added + 1;
            position_in_elements_idx = num_elements+num_elements_added;
            hash_table.elements_idx(position_in_elements_idx) = idx; % Extend element list
            i1 = (position_in_elements_idx-1)*element_DOF+1; % Get the start and end global row indices of the new element
            i2 = (position_in_elements_idx)*element_DOF;
            assert(i2-i1==element_DOF-1);
            fval(i1:i2) = new_element_value; % Extend coefficient list with near zero magnitude (ideally would be zero)
            if refine_previous
                fval_previous(i1:i2) = new_element_value; % Extend coefficient list of previous time step also
            end
            
            [lev_vec, pos_vec] = md_idx_to_lev_pos(num_dims, pde.max_lev, idx);
            
            hash_table.elements.lev_p1(idx,:) = lev_vec+1; % NOTE : have to start lev  index from 1 for sparse storage
            hash_table.elements.pos_p1(idx,:) = pos_vec+1; % NOTE : have to start cell index from 1 for sparse storage
            
        end
        
        %%
        % Set element to type to leaf
        
        [num_live_children, has_complete_children] = ...
            number_of_live_children_idx (hash_table, idx, pde.max_lev, refinement_method);

        if has_complete_children
            hash_table.elements.type(idx) = 1;           
        else
            hash_table.elements.type(idx) = 2;
        end
        
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
if refine_previous
    assert(numel(fval_previous)==numel(hash_table.elements_idx)*element_DOF);
end
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));

for i=1:numel(hash_table.elements_idx)
    lev_vec = hash_table.elements.lev_p1(hash_table.elements_idx(i));
    assert(min(lev_vec>0));
end

%%
% Plot the refined grid (1D only)

plot_grid = 1;
if plot_grid && ~opts.quiet
    coordinates = plot_adapt(pde,opts,hash_table,3);
    plot_adapt_triangle(pde,opts,hash_table,9);
end

elements_idx0 = hash_table.elements_idx;

%%
% Update all the setup outputs which need updating on the new element list

%%
% Update the FMWT transform matrices

for d=1:num_dims
    pde.dimensions{d}.lev = max(hash_table.elements.lev_p1(:,d)-1);
    pde.dimensions{d}.FMWT = OperatorTwoScale(deg,pde.dimensions{d}.lev);
end

%%
% Re check the PDE

pde = check_pde(pde);

%%
% Update the coeff mats to the new size

t = 0;
TD = 0;
pde = get_coeff_mats(pde,t,TD);

%%
% Update A_data

A_data = global_matrix(pde,opts,hash_table);

%%
% Update the conversion to realspace matrices

for d=1:num_dims
    [Meval{d},nodes{d}] = matrix_plot_D(pde,pde.dimensions{d});
end

%%
% Update the coordinates for realspace evaluation

coord = get_realspace_coords(pde,nodes);

%%
% Get the new real space solution and check against unrefined solution

fval_realspace_refined = wavelet_to_realspace(pde,opts,Meval,fval,hash_table);

if ~opts.quiet
    if num_dims == 1
        
        subplot(2,3,4)
        plot(fval)
        hold on
        plot(fval0)
        hold off
        subplot(2,3,5)
        plot(nodes0{1},fval_realspace0)
        hold on
        plot(nodes{1},fval_realspace_refined)
        plot(coordinates,coordinates*0,'o');
        hold off
        
    elseif num_dims == 2
        
        subplot(3,3,4)
        deg1=pde0.deg;
        lev1=pde0.dimensions{1}.lev;
        deg2=pde0.deg;
        lev2=pde0.dimensions{2}.lev;
        dof1=deg1*2^lev1;
        dof2=deg2*2^lev2;
        dofD = dof1*dof2;
        assert(dofD==numel(fval_realspace0));
        f2d = reshape(fval_realspace0,dof2,dof1);
        x = nodes0{1};
        y = nodes0{2};
        if norm(f2d-f2d(1,1))>0 % catch for zero
            contourf(x,y,f2d,'LineColor','none');
        end
        
        subplot(3,3,5)
        deg1=pde.deg;
        lev1=pde.dimensions{1}.lev;
        deg2=pde.deg;
        lev2=pde.dimensions{2}.lev;
        dof1=deg1*2^lev1;
        dof2=deg2*2^lev2;
        dofD = dof1*dof2;
        assert(dofD==numel(fval_realspace_refined));
        f2d = reshape(fval_realspace_refined,dof2,dof1);
        x = nodes{1};
        y = nodes{2};
        if norm(f2d-f2d(1,1))>0 % catch for zero
            contourf(x,y,f2d,'LineColor','none');
        end
        hold on
        scatter(coordinates(:,1),coordinates(:,2),'+','MarkerEdgeColor','white')
        hold off
        
        subplot(3,3,6)
        fval_element = zeros(num_elements,1);
        for i=1:deg^num_dims:num_elements
            view = fval((i-1).*element_DOF+1:i*element_DOF);
            fval_element(i) = sqrt(sum(view.^2));
        end
        tmp = nonzeros(fval_element);
        semilogy(tmp);
        
    end
end

assert(numel(fval)==numel(hash_table.elements_idx)*element_DOF);
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));
assert(sum(hash_table.elements_idx-elements_idx0)==0);

end
