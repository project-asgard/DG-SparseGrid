function [pde,fval,hash_table,A_data,Meval,nodes,coord] = adapt(pde,opts,fval,hash_table,nodes0,fval_realspace0)

num_elements    = numel(hash_table.elements_idx);
num_dimensions  = numel(pde.dimensions);
deg             = pde.deg; 

element_DOF = deg^num_dimensions;

refine_threshold  = max(abs(fval)) * 1e-2;
coarsen_threshold = max(abs(fval)) * 1e-4;

newElementVal = 1e-15;

debug   = 0;
coarsen = 1;
refine  = 1;
method  = 1;

assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));

if ~opts.quiet
    fprintf('Initial number of elementss: %i\n', numel(hash_table.elements_idx));
    fprintf('Initial number of DOFs: %i\n', numel(hash_table.elements_idx)*deg^num_dimensions);
end

%%
% Store unrefined fval for comparison with refined fval after

pde0  = pde;
fval0 = fval;

%%
% Plot the grid (1D only)

plot_grid = 1;
if plot_grid && ~opts.quiet
    plot_adapt(pde,opts,hash_table,1);
    plot_adapt_triangle(pde,opts,hash_table,7);
end

%%
% Coarsen

if coarsen
    
    num_remove = 0;
    remove_list = [];
    
    for n=1:num_elements
        
        idx = hash_table.elements_idx(n);
        
        gidx1 = (n-1)*element_DOF+1;
        gidx2 = gidx1 + element_DOF - 1;
        
        element_sum = sum(abs(fval(gidx1:gidx2)),'all');
                
        %%
        % Check for coarsening (de-refinement)
        
        if element_sum <= coarsen_threshold % Check only the deg=0 term for each element
            
            if debug; fprintf('leaf node to be REMOVED, fval=%f\n', element_sum); end
            
            thisElemLevVec = hash_table.elements.lev_p1(idx,:)-1; % NOTE : remove the 1 per note below
            thisElemPosVec = hash_table.elements.pos_p1(idx,:)-1; % NOTE : remove the 1 per note below
            
            %%
            % Generate a list of elements who will become leaves after
            % removing this element.
            
            for d=1:num_dimensions
                
                newLeafElemLevVec = thisElemLevVec;
                newLeafElemPosVec = thisElemPosVec;
                
                newLeafElemLevVec(d) = newLeafElemLevVec(d)-1;
                newLeafElemPosVec(d) = floor(newLeafElemPosVec(d)/2);
                
                idx = lev_cell_to_element_index(pde,newLeafElemLevVec,newLeafElemPosVec);
                
                %%
                % Assert this element exists
                
                %                     assert(hash_table.elements.lev_p1(element_idx,d) ~= 0);
                
                %%
                % Set element type to leaft == 2
                
                %                     hash_table.elements.node_type(element_idx) = 2;
                
            end
            
            num_remove = num_remove + 1;
            remove_list(num_remove) = n;
            
        end
    end
    
    %%
    % Now remove elements
    
    assert(numel(remove_list)==num_remove);
    
    remove_list2 = [];
    for n=1:num_remove
        
        %%
        % Remove entries from element table (recall sparse storage means =0
        % removes it from the table
        
        hash_table.elements.lev_p1(hash_table.elements_idx(remove_list(n)),:) = 0;
        hash_table.elements.pos_p1(hash_table.elements_idx(remove_list(n)),:) = 0;
        hash_table.elements.node_type(hash_table.elements_idx(remove_list(n)))= 0;
        
        %%
        % Remove all deg parts of this element from fval
        
        nn = remove_list(n);
        
        i1 = (nn-1)*element_DOF+1; % Get the start and end global row indices of the element
        i2 = (nn)*element_DOF;
        assert(i2-i1==element_DOF-1);
        
        remove_list2 = [remove_list2,i1:i2];
        
    end
    
    %%
    % Remove entries from elements_idx
    
    hash_table.elements_idx(remove_list) = [];
    fval(remove_list2) = [];
    if ~opts.quiet; fprintf('    Coarsen on : removed %i elements\n', num_remove); end

end

assert(numel(fval)==numel(hash_table.elements_idx)*element_DOF);
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));

%%
% Plot the refined grid (1D only)

plot_grid = 1;
if plot_grid && ~opts.quiet
    plot_adapt(pde,opts,hash_table,2);
    plot_adapt_triangle(pde,opts,hash_table,8);
end

%%
% Refine

if refine
    
    num_elements = numel(hash_table.elements_idx);
    cnt = 0;
    clear newElemLevVecs;
    clear newElemPosVecs;
    for n=1:num_elements
        
        idx = hash_table.elements_idx(n);
        lev_vec = hash_table.elements.lev_p1(idx,:)-1;
        pos_vec = hash_table.elements.pos_p1(idx,:)-1;
        
        gidx1 = (n-1)*element_DOF+1;
        gidx2 = gidx1 + element_DOF - 1;
        
        element_sum = sum(abs(fval(gidx1:gidx2)),'all');
               
        %%
        % Check for refinement
        
        if element_sum >= refine_threshold
            
            if debug; fprintf('leaf node to be refined, fval=%f\n', element_sum); end
            
            [daughterElemLevVecs,daughterElemPosVecs,nDaughters] = get_my_daughters(lev_vec, pos_vec, pde.max_lev, method);
            
            newElemLevVecs(cnt+1:cnt+nDaughters,:) = daughterElemLevVecs;
            newElemPosVecs(cnt+1:cnt+nDaughters,:) = daughterElemPosVecs;
            
            if nDaughters > 0
                hash_table.elements.node_type(idx) = 1; % Now that this element has been refined it is no longer a leaf.              
            end
            
            cnt = cnt + nDaughters;
            
        else
            
            if debug; fprintf('leaf node but no refinement, fval=%f\n', element_sum); end
            
        end
        
    end
    
    %%
    % Now add these elements with (almost) zero coefficient to the
    % elements table and elementsIDX
    
    num_try_to_add = cnt;
    num_add = 0;
    for i=1:num_try_to_add
        
        thisElemLevVec = newElemLevVecs(i,:);
        thisElemPosVec = newElemPosVecs(i,:);
        idx = lev_cell_to_element_index(thisElemLevVec,thisElemPosVec,pde.max_lev);
        assert(idx>=0);
        assert(min(thisElemLevVec)>=0);
        assert(min(thisElemPosVec)>=0);
        
        %%
        % If element does not exist, add it
        
        if hash_table.elements.lev_p1(idx,1)==0
            
            num_add = num_add + 1;
            myIdx = num_elements+num_add;
            hash_table.elements_idx(myIdx) = idx; % Extend element list
            i1 = (myIdx-1)*element_DOF+1; % Get the start and end global row indices of the new element
            i2 = (myIdx)*element_DOF;
            assert(i2-i1==element_DOF-1);
            fval(i1:i2) = newElementVal; % Extend coefficient list with near zero magnitude (ideally would be zero)
            
            hash_table.elements.lev_p1(idx,:) = thisElemLevVec+1; % NOTE : have to start lev  index from 1 for sparse storage
            hash_table.elements.pos_p1(idx,:) = thisElemPosVec+1; % NOTE : have to start cell index from 1 for sparse storage
            
        end
        
        %%
        % Set element to type to leaf
        
        hash_table.elements.node_type(idx) = 2;
        
    end
    
    if ~opts.quiet; fprintf('    Refine on : added %i elements\n', num_add); end
    
end
assert(numel(fval)==numel(hash_table.elements_idx)*element_DOF);
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));

for i=1:numel(hash_table.elements_idx)
    lev_vec = hash_table.elements.lev_p1(hash_table.elements_idx(i));
    assert(min(lev_vec>0));
end

%%
% Plot the refined grid (1D only)

plot_grid = 1;
if plot_grid && ~opts.quiet
    plot_adapt(pde,opts,hash_table,3);
    plot_adapt_triangle(pde,opts,hash_table,9);
end

elements_idx0 = hash_table.elements_idx;

%%
% Update all the setup outputs which need updating on the new element list

%%
% Update the FMWT transform matrices

for d=1:num_dimensions
    pde.dimensions{d}.lev = max(hash_table.elements.lev_p1(:,d)-1);
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde,d,deg,pde.dimensions{d}.lev);
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

for d=1:num_dimensions
    [Meval{d},nodes{d}] = matrix_plot_D(pde,pde.dimensions{d});
end

%%
% Update the coordinates for realspace evaluation

coord = get_realspace_coords(pde,nodes);

%%
% Get the new real space solution and check against unrefined solution

fval_realspace_refined = multi_2D_D(pde,opts,Meval,fval,hash_table);

if ~opts.quiet
    if num_dimensions == 1
        
        subplot(2,3,4)
        plot(fval)
        hold on
        plot(fval0)
        hold off
        subplot(2,3,5)
        plot(nodes0{1},fval_realspace0)
        hold on
        plot(nodes{1},fval_realspace_refined)
        hold off
        
    elseif num_dimensions == 2
        
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
        contour(x,y,f2d);
        
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
        contour(x,y,f2d);
        
    end
end

assert(numel(fval)==numel(hash_table.elements_idx)*element_DOF);
assert(numel(find(hash_table.elements_idx))==numel(hash_table.elements_idx));
assert(sum(hash_table.elements_idx-elements_idx0)==0);

if ~opts.quiet
    fprintf('Final number of elementss: %i\n', numel(hash_table.elements_idx));
    fprintf('Final number of DOFs: %i\n', numel(hash_table.elements_idx)*deg^num_dimensions);
end

end