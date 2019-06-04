function [pde,fval,A_data,Meval,nodes,coord] = adapt(pde,opts,fval,HASHInv,connectivity,nodes0,fval_realspace0)

num_elements    = numel(pde.elementsIDX);
num_dimensions  = numel(pde.dimensions);
deg             = pde.dimensions{1}.deg; % TODO

elementDOF = deg^num_dimensions;

refine_threshold  = max(abs(fval)) * 1e-4;
coarsen_threshold = max(abs(fval)) * 1e-6;

newElementVal = 1e-15;

debug   = 0;
coarsen = 1;
refine  = 1;

assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));

fprintf('Initial number of elementss: %i\n', numel(pde.elementsIDX));
fprintf('Initial number of DOFs: %i\n', numel(pde.elementsIDX)*deg^num_dimensions);

%%
% Store unrefined fval for comparison with refined fval after

pde0  = pde;
fval0 = fval;

%%
% Plot the grid (1D only)

plot_grid = 1;
if plot_grid
    plot_adapt(pde,1);
    plot_adapt_triangle(pde,7);
end

%%
% Coarsen

if coarsen
    
    num_remove = 0;
    remove_list = [];
    
    for n=1:num_elements
        
        idx = pde.elementsIDX(n);
        
        gidx1 = (n-1)*elementDOF+1;
        gidx2 = gidx1 + elementDOF - 1;
        
        element_sum = sum(abs(fval(gidx1:gidx2)),'all');
        
%         if pde.elements.node_type(idx) == 2 % refine leaf nodes
            
            %%
            % Check for coarsening (de-refinement)
            
            if element_sum <= coarsen_threshold % Check only the deg=0 term for each element
                
                if debug; fprintf('leaf node to be REMOVED, fval=%f\n', element_sum); end
                
                thisElemLevVec = pde.elements.lev_p1(idx,:)-1; % NOTE : remove the 1 per note below
                thisElemPosVec = pde.elements.pos_p1(idx,:)-1; % NOTE : remove the 1 per note below
                
                %%
                % Generate a list of elements who will become leaves after
                % removing this element.
                
                for d=1:num_dimensions
                    
                    newLeafElemLevVec = thisElemLevVec;
                    newLeafElemPosVec = thisElemPosVec;
                    
                    newLeafElemLevVec(d) = newLeafElemLevVec(d)-1;
                    newLeafElemPosVec(d) = floor(newLeafElemPosVec(d)/2);
                    
                    element_idx = lev_cell_to_element_index(pde,newLeafElemLevVec,newLeafElemPosVec);
                    
                    %%
                    % Assert this element exists
                    
%                     assert(pde.elements.lev_p1(element_idx,d) ~= 0);
                    
                    %%
                    % Set element type to leaft == 2
                    
%                     pde.elements.node_type(element_idx) = 2;
                    
                end
                
                num_remove = num_remove + 1;
                remove_list(num_remove) = n;
                
            end
        end
%     end
    
    %%
    % Now remove elements
    
    assert(numel(remove_list)==num_remove);
    
    remove_list2 = [];
    for n=1:num_remove
        
        %%
        % Remove entries from element table (recall sparse storage means =0
        % removes it from the table
        
        pde.elements.lev_p1(pde.elementsIDX(remove_list(n)),:) = 0;
        pde.elements.pos_p1(pde.elementsIDX(remove_list(n)),:) = 0;
        pde.elements.node_type(pde.elementsIDX(remove_list(n)))= 0;
        
        %%
        % Remove all deg parts of this element from fval
        
        nn = remove_list(n);
        
        i1 = (nn-1)*elementDOF+1; % Get the start and end global row indices of the element
        i2 = (nn)*elementDOF;
        assert(i2-i1==elementDOF-1);
        
        remove_list2 = [remove_list2,i1:i2];
        
    end
    
    %%
    % Remove entries from elementsIDX
    
    pde.elementsIDX(remove_list) = [];
    fval(remove_list2) = [];
    
end

fprintf('    Removed %i elements\n', num_remove);

assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);
assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));

%%
% Plot the refined grid (1D only)

plot_grid = 1;
if plot_grid
   plot_adapt(pde,2);
   plot_adapt_triangle(pde,8);
end

%%
% Refine

if refine
    
    num_elements = numel(pde.elementsIDX);
    cnt = 0;
    clear newElemLevVecs;
    clear newElemPosVecs;
    for n=1:num_elements
        
        idx = pde.elementsIDX(n);
        lev_vec = pde.elements.lev_p1(idx,:)-1;
        pos_vec = pde.elements.pos_p1(idx,:)-1;
        
        gidx1 = (n-1)*elementDOF+1;
        gidx2 = gidx1 + elementDOF - 1;
        
        element_sum = sum(abs(fval(gidx1:gidx2)),'all');
        
%         if pde.elements.node_type(idx) == 2 % refine leaf nodes only according to their deg=0 element
            
            %%
            % Check for refinement
            
%             fprintf('element_sum: %f\n', element_sum);
            
            if element_sum >= refine_threshold
                
                if debug; fprintf('leaf node to be refined, fval=%f\n', element_sum); end
                                
                [daughterElemLevVecs,daughterElemPosVecs,nDaughters] = get_my_daughters(lev_vec,pos_vec,pde.maxLev);
                
                newElemLevVecs(cnt+1:cnt+nDaughters,:) = daughterElemLevVecs;
                newElemPosVecs(cnt+1:cnt+nDaughters,:) = daughterElemPosVecs; 

                if nDaughters > 0
                    pde.elements.node_type(idx) = 1; % Now that this element has been refined it is no longer a leaf.

                end
                
                cnt = cnt + nDaughters;
                                
            else
                
                if debug; fprintf('leaf node but no refinement, fval=%f\n', element_sum); end
                
            end
            
%         else
            
            if debug; disp('internal node'); end
            
%         end
        
    end
    
    %%
    % Now add these elements with (almost) zero coefficient to the
    % elements table and elementsIDX
    
    num_add = cnt;
    add_cnt = 0;
    for i=1:num_add
        
        thisElemLevVec = newElemLevVecs(i,:);
        thisElemPosVec = newElemPosVecs(i,:);
        element_idx = lev_cell_to_element_index(pde,thisElemLevVec,thisElemPosVec);
        assert(element_idx>=0);
        
        %%
        % If element does not exist, add it
        
        if pde.elements.lev_p1(element_idx,1)==0
            
            add_cnt = add_cnt + 1;
            myIdx = num_elements+add_cnt;
            pde.elementsIDX(myIdx) = element_idx; % Extend element list
            i1 = (myIdx-1)*elementDOF+1; % Get the start and end global row indices of the new element
            i2 = (myIdx)*elementDOF;
            assert(i2-i1==elementDOF-1);
            fval(i1:i2) = newElementVal; % Extend coefficient list with near zero magnitude (ideally would be zero)
            
            pde.elements.lev_p1(element_idx,:) = thisElemLevVec+1; % NOTE : have to start lev  index from 1 for sparse storage
            pde.elements.pos_p1(element_idx,:) = thisElemPosVec+1; % NOTE : have to start cell index from 1 for sparse storage
            
        end
        
        %%
        % Set element to type to leaf
        
        pde.elements.node_type(element_idx) = 2;
        
    end
end
assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);
assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));

fprintf('    Added %i elements\n', add_cnt);


%%
% Set any leafs with all leaf daughters to be a node

leafCheck = 0;
if leafCheck
    for n=1:numel(pde.elementsIDX)
        if pde.elements.node_type(pde.elementsIDX(n)) == 2
            idx = pde.elementsIDX(n);
            lev_vec = pde.elements.lev_p1(idx,:)-1;
            pos_vec = pde.elements.pos_p1(idx,:)-1;
            [daughterElemLevVecs,daughterElemPosVecs,nDaughters] = get_my_daughters(lev_vec,pos_vec,pde.maxLev);
            
            nLeaves = 0;
            for nn=1:nDaughters
                thisElemLevVec = daughterElemLevVecs(nn,:);
                thisElemPosVec = daughterElemPosVecs(nn,:);
                idx2 = lev_cell_to_element_index(pde,thisElemLevVec,thisElemPosVec);
                if pde.elements.node_type(idx2) ~= 0 % check for existence
                    nLeaves = nLeaves + 1;
                end
            end
            
            if nLeaves == nDaughters && nDaughters > 0
                pde.elements.node_type(idx) = 1; % No longer a leaf
            end
            
        end
    end
end
assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);
assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));


%%
% Plot the refined grid (1D only)

plot_grid = 1;
if plot_grid
   plot_adapt(pde,3);
   plot_adapt_triangle(pde,9);
end

elementsIDX0 = pde.elementsIDX;

%%
% Update all the setup outputs which need updating on the new element list

%%
% Update the FMWT transform matrices

for d=1:num_dimensions
    pde.dimensions{d}.lev = max(pde.elements.lev_p1(:,d)-1);
    pde.dimensions{d}.FMWT = OperatorTwoScale(pde,d,deg,pde.dimensions{d}.lev);
end

%%
% Update the coeff mats to the new size

t = 0;
TD = 0;
pde = getCoeffMats(pde,t,TD);

%%
% Update A_data

A_data = GlobalMatrixSG_SlowVersion(pde,opts,HASHInv,connectivity,deg);

%%
% Update the conversion to realspace matrices

for d=1:num_dimensions
    [Meval{d},nodes{d}] = matrix_plot_D(pde.dimensions{d});
end

%%
% Update the coordinates for realspace evaluation

coord = get_realspace_coords(pde,nodes);

%%
% Get the new real space solution and check against unrefined solution

fval_realspace_refined = Multi_2D_D(pde,Meval,fval,HASHInv);

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
    deg1=pde0.dimensions{1}.deg;
    lev1=pde0.dimensions{1}.lev;
    deg2=pde0.dimensions{2}.deg;
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
    deg1=pde.dimensions{1}.deg;
    lev1=pde.dimensions{1}.lev;
    deg2=pde.dimensions{2}.deg;
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

assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);
assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));
assert(sum(pde.elementsIDX-elementsIDX0)==0);

fprintf('Final number of elementss: %i\n', numel(pde.elementsIDX));
fprintf('Final number of DOFs: %i\n', numel(pde.elementsIDX)*deg^num_dimensions);

end