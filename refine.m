function [pde,fval,A_data,Meval,nodes,coord] = refine(pde,opts,fval,HASHInv,connectivity)

nElements = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);
deg = pde.dimensions{1}.deg; % TODO

elementDOF = deg^nDims;

threshold = 1e-5;
rel_threshold = max(fval)*threshold;
newElementVal = 1e-15;

debug = 0;

%%
% Plot the grid (1D only)

plot_grid = 1;
if plot_grid
    if nDims == 1
        figure(222);
        hold on
        for i=1:nElements
            x = pde.elements.pos(pde.elementsIDX(i))-1;
            y = pde.elements.lev(pde.elementsIDX(i))-1;
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

cnt = 1;
for n=1:nElements
    
    idx = pde.elementsIDX(n);
    gidx = (n-1)*elementDOF+1;
    
    if pde.elements.node_type(idx) == 2 % refine leaf nodes
        
        if abs(fval(gidx)) >= rel_threshold
            
            if debug; fprintf('leaf node to be refined, fval=%f\n',fval(gidx)); end
            
            %%
            % Get this coordinate vector
            % levVec = [lev1,lev2,...,levD]
            % posVec = [pos1,pos2,...,posD]
            
            thisElemLevVec = pde.elements.lev(idx,:)-1; % NOTE : remove the 1 per note below
            thisElemPosVec = pde.elements.pos(idx,:)-1; % NOTE : remove the 1 per note below
            
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
                    
                    newElemLevVecs(cnt,:) = thisElemLevVec;
                    newElemPosVecs(cnt,:) = thisElemPosVec;
                    
                    newElemLevVecs(cnt,d2) = thisElemLevVec(d)+1;
                    newElemPosVecs(cnt,d2) = thisElemPosVec(d)*2; % Assumes pos starts at 0
                                   
                    assert(newElemPosVecs(cnt,d2) >= 0);
                    assert(newElemLevVecs(cnt,d2) >= 0);
                                     
                    cnt = cnt + 1;
                    
                    %%
                    % Second daughter
                    
                    newElemLevVecs(cnt,:) = thisElemLevVec;
                    newElemPosVecs(cnt,:) = thisElemPosVec;
                    
                    newElemLevVecs(cnt,d2) = thisElemLevVec(d)+1;
                    newElemPosVecs(cnt,d2) = thisElemPosVec(d)*2+1; % Assumes pos starts at 0
                    
                    cnt = cnt + 1;
                                        
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
                pde.elementsIDX(nElements+i) = element_idx; % Extend element list
                i1 = (nElements+i-1)*elementDOF-1; % Get the start and end global row indices of the new element
                i2 = (nElements+i)*elementDOF;
                fval(i1:i2) = newElementVal; % Extend coefficient list with near zero magnitude (ideally would be zero)
                
                pde.elements.lev(element_idx,:) = thisElemLevVec+1; % NOTE : have to start lev  index from 1 for sparse storage
                pde.elements.pos(element_idx,:) = thisElemPosVec+1; % NOTE : have to start cell index from 1 for sparse storage
                pde.elements.node_type(element_idx) = 2; % This new element is now a leaf
                
            end
            
            pde.elements.node_type(idx) = 1; % Now that this element has been refined it is no longer a leaf.
            
        else
            
            if debug; fprintf('leaf node but no refinement, fval=%f\n',fval(gidx)); end
        end
        
    else
        
        if debug; disp('internal node'); end
        
    end
    
end

assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);

%%
% Update all the setup outputs which need updating on the new element list

%%
% Update the FMWT transform matrices

for d=1:nDims
    pde.dimensions{d}.lev = max(pde.elements.lev(:,d)-1);
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