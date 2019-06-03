function [pde,fval,A_data,Meval,nodes,coord] = adapt(pde,opts,fval,HASHInv,connectivity,nodes0,fval_realspace0)

num_elements = numel(pde.elementsIDX);
num_dims     = numel(pde.dimensions);
deg          = pde.dimensions{1}.deg; % TODO

elementDOF = deg^num_dims;

refine_threshold  = max(abs(fval)) * 1e-6;
coarsen_threshold = max(abs(fval)) * 1e-10;

newElementVal = 1e-15;

debug   = 1;
coarsen = 1;
refine  = 1;

assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));

%%
% Store unrefined fval for comparison with refined fval after

pde0  = pde;
fval0 = fval;

%%
% Plot the grid (1D only)

plot_grid = 1;
if plot_grid
    plot_adapt(pde,1);
end

%%
% Coarsen

if coarsen
    
    nRemove = 0;
    remove_list = [];
    
    for n=1:num_elements
        
        idx = pde.elementsIDX(n);
%         gidx = (n-1)*elementDOF+1; % this is the deg=0 part of this element
        
        gidx1 = (n-1)*elementDOF+1;
        gidx2 = gidx1 + elementDOF - 1;
        
        element_sum = sum(abs(fval(gidx1:gidx2)),'all');
        
        if pde.elements.node_type(idx) == 2 % refine leaf nodes
            
            %%
            % Check for coarsening (de-refinement)
            
            if element_sum <= coarsen_threshold % Check only the deg=0 term for each element
                
                if debug; fprintf('leaf node to be REMOVED, fval=%f\n', element_sum); end
                
                thisElemLevVec = pde.elements.lev_p1(idx,:)-1; % NOTE : remove the 1 per note below
                thisElemPosVec = pde.elements.pos_p1(idx,:)-1; % NOTE : remove the 1 per note below
                
                %%
                % Generate a list of elements who will become leaves after
                % removing this element.
                
                for d=1:num_dims
                    
                    newLeafElemLevVec = thisElemLevVec;
                    newLeafElemPosVec = thisElemPosVec;
                    
                    newLeafElemLevVec(d) = newLeafElemLevVec(d)-1;
                    newLeafElemPosVec(d) = floor(newLeafElemPosVec(d)/2);
                    
                    element_idx = lev_cell_to_element_index(pde,newLeafElemLevVec,newLeafElemPosVec);
                    
                    %%
                    % Assert this element exists
                    
                    assert(pde.elements.lev_p1(element_idx,d) ~= 0);
                    
                    %%
                    % Set element type to leaft == 2
                    
                    pde.elements.node_type(element_idx) = 2;
                    
                end
                
                nRemove = nRemove + 1;
                remove_list(nRemove) = n;
                
            end
        end
    end
    
    %%
    % Now remove elements
    
    assert(numel(remove_list)==nRemove);
    
    remove_list2 = [];
    for n=1:nRemove
        
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
assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);
assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));

%%
% Plot the refined grid (1D only)

plot_grid = 1;
if plot_grid
   plot_adapt(pde,2);
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
%         gidx = (n-1)*elementDOF+1; % this is deg=0 part of this element
        
        gidx1 = (n-1)*elementDOF+1;
        gidx2 = gidx1 + elementDOF - 1;
        
        element_sum = sum(abs(fval(gidx1:gidx2)),'all');
        
        if pde.elements.node_type(idx) == 2 % refine leaf nodes only according to their deg=0 element
            
            %%
            % Check for refinement
            
            if element_sum >= refine_threshold
                
                if debug; fprintf('leaf node to be refined, fval=%f\n', element_sum); end
                
                nDaughters = 0;
                [daughterElemLevVecs,daughterElemPosVecs,nDaughters] = getMyDaughters(pde,idx);
                
                newElemLevVecs(cnt+1:cnt+nDaughters,:) = daughterElemLevVecs;
                newElemPosVecs(cnt+1:cnt+nDaughters,:) = daughterElemPosVecs; 

                if nDaughters > 0
                    pde.elements.node_type(idx) = 1; % Now that this element has been refined it is no longer a leaf.

                end
                
                cnt = cnt + nDaughters;
                                
            else
                
                if debug; fprintf('leaf node but no refinement, fval=%f\n',fval(gidx)); end
                
            end
            
        else
            
            if debug; disp('internal node'); end
            
        end
        
    end
    
    %%
    % Now add these elements with (almost) zero coefficient to the
    % elements table and elementsIDX
    
    nAdd = cnt;
    addCnt = 0;
    for i=1:nAdd
        
        thisElemLevVec = newElemLevVecs(i,:);
        thisElemPosVec = newElemPosVecs(i,:);
        element_idx = lev_cell_to_element_index(pde,thisElemLevVec,thisElemPosVec);
        assert(element_idx>=0);
        
        %%
        % If element does not exist, add it
        
        if pde.elements.lev_p1(element_idx,1)==0
            
            addCnt = addCnt + 1;
            myIdx = num_elements+addCnt;
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
        
        if thisElemLevVec(1) == 3 && thisElemPosVec(1)==2
            disp('');
        end
        
    end
end
assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);
assert(numel(find(pde.elementsIDX))==numel(pde.elementsIDX));

%%
% Set any leafs with all leaf daughters to be a node

leafCheck = 1;
if leafCheck
    for n=1:numel(pde.elementsIDX)
        if pde.elements.node_type(pde.elementsIDX(n)) == 2
            idx = pde.elementsIDX(n);
            [daughterElemLevVecs,daughterElemPosVecs,nDaughters] = getMyDaughters(pde,idx);
            
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
end

elementsIDX0 = pde.elementsIDX;

%%
% Update all the setup outputs which need updating on the new element list

%%
% Update the FMWT transform matrices

for d=1:num_dims
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

for d=1:num_dims
    [Meval{d},nodes{d}] = matrix_plot_D(pde.dimensions{d});
end

%%
% Update the coordinates for realspace evaluation

coord = get_realspace_coords(pde,nodes);

%%
% Get the new real space solution and check against unrefined solution

fval_realspace_refined = Multi_2D_D(pde,Meval,fval,HASHInv);

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
    hold off
elseif num_dims == 2
    
    subplot(2,3,4)
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
    
    subplot(2,3,5)
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

end