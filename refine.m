function [pde,fval,A_data,Meval,nodes,coord] = refine(pde,opts,fval,HASHInv,connectivity)

nElements = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);
deg = pde.dimensions{1}.deg; % TODO

elementDOF = deg^nDims;

refine_threshold  = max(fval) * 1e-5;
coarsen_threshold = max(fval) * 1e-7;

newElementVal = 1e-15;

debug = 1;
coarsen = 1;
refine = 0;

%%
% Plot the grid (1D only)

plot_grid = 1;
if plot_grid
    if nDims == 1
        figure(222);
        subplot(1,2,1)
        hold on
        for i=1:nElements
            x = pde.elements.pos_p1(pde.elementsIDX(i))-1;
            y = pde.elements.lev_p1(pde.elementsIDX(i))-1;
            c = pde.elements.node_type(pde.elementsIDX(i));
            if c == 1
                style = 'ob';
            elseif c == 2
                style = 'or';
            end
            offset = 2^(y-1)/2;
            if y > 1
                
                s = 2^(y-1)-1;
                h = 1/(2^(y-1));
                w = 1-h;
                o = h/2;
                plot(x/s*w+o,-y,style);
                
            else
                plot(x+0.5,-y,style);
            end
        end
        hold off
    end
end

%%
% Determine which elements to refine and add them.

if coarsen
    
    nRemove = 0;
    remove_list = [];
    
    for n=1:nElements
        
        idx = pde.elementsIDX(n);
        gidx = (n-1)*elementDOF+1;
        
        if pde.elements.node_type(idx) == 2 % refine leaf nodes
            
            %%
            % Check for coarsening (de-refinement)
            
            if abs(fval(gidx)) <= coarsen_threshold
                
                if debug; fprintf('leaf node to be REMOVED, fval=%f\n',fval(gidx)); end
                
                thisElemLevVec = pde.elements.lev_p1(idx,:)-1; % NOTE : remove the 1 per note below
                thisElemPosVec = pde.elements.pos_p1(idx,:)-1; % NOTE : remove the 1 per note below
                
                %%
                % Generate a list of elements who will become leaves after
                % removing this element.
                
                for d=1:nDims
                    
                    newLeafElemLevVec = thisElemLevVec;
                    newLeafElemPosVec = thisElemPosVec;
                    
                    newLeafElemLevVec(d) = newLeafElemLevVec(d)-1;
                    newLeafElemPosVec(d) = floor(newLeafElemPosVec(d)/2);
                    
                    element_idx = lev_cell_to_element_index(pde,newLeafElemLevVec,newLeafElemPosVec);
                    
                    %%
                    % Assert this element exists
                    
                    assert(pde.elements.lev_p1(element_idx,1) ~= 0);
                    
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
    for n=1:nRemove
        
        %%
        % Remove entries from element table
        
        pde.elements.lev_p1(pde.elementsIDX(remove_list(n)),:) = 0;
        pde.elements.pos_p1(pde.elementsIDX(remove_list(n)),:) = 0;
        pde.elements.node_type(pde.elementsIDX(remove_list(n)))= 0;
        
        %%
        % Remove enetries from fval
        
        i1 = (n-1)*elementDOF+1; % Get the start and end global row indices of the new element
        i2 = (n)*elementDOF;
        assert(i2-i1==elementDOF-1);
        
        fval(i1:i2) = []; % Extend coefficient list with near zero magnitude (ideally would be zero)
        
    end
    
    %%
    % Remove entries from elementsIDX
    
    pde.elementsIDX(remove_list) = [];
    
end
assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);


if refine
    
    nElements = numel(pde.elementsIDX);
    cnt = 1;
    clear newElemLevVecs;
    clear newElemPosVecs;
    for n=1:nElements
        
        idx = pde.elementsIDX(n);
        gidx = (n-1)*elementDOF+1;
        
        if pde.elements.node_type(idx) == 2 % refine leaf nodes
            
            %%
            % Check for refinement
            
            if abs(fval(gidx)) >= refine_threshold
                
                if debug; fprintf('leaf node to be refined, fval=%f\n',fval(gidx)); end
                
                %%
                % Get this coordinate vector
                % levVec = [lev1,lev2,...,levD]
                % posVec = [pos1,pos2,...,posD]
                
                thisElemLevVec = pde.elements.lev_p1(idx,:)-1; % NOTE : remove the 1 per note below
                thisElemPosVec = pde.elements.pos_p1(idx,:)-1; % NOTE : remove the 1 per note below
                
                %%
                % Generate list of new elements which satisfy selection rule
                
                for d=1:nDims
                    
                    %%
                    % Add two daughter nodes in each dimension for each leaf
                    % element
                    
                    for d2=1:nDims
                        
                        %%
                        % TODO : This refinement is not sparse. Need to add the
                        % sparse grid selection rule here also.
                        
                        %%
                        % First daughter
                        
                        newElemLevVec = thisElemLevVec;
                        newElemPosVec = thisElemPosVec;
                        newElemLevVec(d2) = newElemLevVec(d)+1;
                        newElemPosVec(d2) = newElemPosVec(d)*2; % Assumes pos starts at 0
                        
                        if sum(newElemLevVec)<=thisElemLevVec(d)+1 && newElemLevVec(d2)<=pde.maxLev % Sparse grid selection rule AND max depth check
                            
                            newElemLevVecs(cnt,:) = newElemLevVec;
                            newElemPosVecs(cnt,:) = newElemPosVec; % Assumes pos starts at 0
                            
                            assert(newElemPosVecs(cnt,d2) >= 0);
                            assert(newElemLevVecs(cnt,d2) >= 0);
                            
                            cnt = cnt + 1;
                            
                            pde.elements.node_type(idx) = 1; % Now that this element has been refined it is no longer a leaf.
                            
                        else
                            
                            disp('element not added because it did not obey sparse selection rule');
                            
                        end
                        
                        %%
                        % Second daughter
                        
                        newElemLevVec = thisElemLevVec;
                        newElemPosVec = thisElemPosVec;
                        
                        newElemLevVec(d2) = newElemLevVec(d)+1;
                        newElemPosVec(d2) = newElemPosVec(d)*2+1; % Assumes pos starts at 0
                        
                        if sum(newElemLevVec)<=thisElemLevVec(d)+1 && newElemLevVec(d2)<=pde.maxLev % Sparse grid selection rule AND max depth check
                            
                            newElemLevVecs(cnt,:) = newElemLevVec;
                            newElemPosVecs(cnt,:) = newElemPosVec; % Assumes pos starts at 0
                            
                            assert(newElemPosVecs(cnt,d2) >= 0);
                            assert(newElemLevVecs(cnt,d2) >= 0);
                            
                            cnt = cnt + 1;
                            
                            pde.elements.node_type(idx) = 1; % Now that this element has been refined it is no longer a leaf.

                        else
                            
                            disp('element not added because it did not obey sparse selection rule');
                            
                            
                        end
                        
                    end
                    
                end
                                
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
    
    nAdd = cnt-1;
    for i=1:nAdd
        
        thisElemLevVec = newElemLevVecs(i,:);
        thisElemPosVec = newElemPosVecs(i,:);
        element_idx = lev_cell_to_element_index(pde,thisElemLevVec,thisElemPosVec);
        
        %%
        % If element does not exist, add it
        
        if pde.elements.lev_p1(element_idx,1)==0
            
            pde.elementsIDX(nElements+i) = element_idx; % Extend element list
            i1 = (nElements+i-1)*elementDOF+1; % Get the start and end global row indices of the new element
            i2 = (nElements+i)*elementDOF;
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

%%
% Plot the refined grid (1D only)

nElements = numel(pde.elementsIDX);
plot_grid = 1;
if plot_grid
    if nDims == 1
        figure(222);
        subplot(1,2,2)
        hold on
        for i=1:nElements
            x = pde.elements.pos_p1(pde.elementsIDX(i))-1;
            y = pde.elements.lev_p1(pde.elementsIDX(i))-1;
            c = pde.elements.node_type(pde.elementsIDX(i));
            if c == 1
                style = 'ob';
            elseif c == 2
                style = 'or';
            end
            offset = 2^(y-1)/2;
            if y > 1
                
                s = 2^(y-1)-1;
                h = 1/(2^(y-1));
                w = 1-h;
                o = h/2;
                plot(x/s*w+o,-y,style);
                
            else
                plot(x+0.5,-y,style);
            end
        end
        hold off
    end
end

%%
% Update all the setup outputs which need updating on the new element list

%%
% Update the FMWT transform matrices

for d=1:nDims
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

for d=1:nDims
    [Meval{d},nodes{d}] = matrix_plot_D(pde.dimensions{d});
end

%%
% Update the coordinates for realspace evaluation

coord = get_realspace_coords(pde,nodes);

end