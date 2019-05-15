function [pde,fval,A_data] = refine(pde,opts,fval,HASHInv,connectivity)

nElements = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);
deg = pde.dimensions{1}.deg; % TODO

elementDOF = deg^nDims;

threshold = 1e-6;
rel_threshold = max(fval)*threshold;
newElementVal = 1e-15;

cnt = 1;
for n=1:nElements
    
    idx = pde.elementsIDX(n);
    gidx = (n-1)*elementDOF+1;
    
    if pde.elements.node_type(idx) == 2 % refine leaf nodes
        
        if abs(fval(gidx)) >= rel_threshold
            
            fprintf('leaf node to be refined, fval=%f\n',fval(gidx))
            
            %%
            % Get this coordinate vector
            % levVec = [lev1,lev2,...,levD]
            % posVec = [pos1,pos2,...,posD]
            
            thisElemLevVec = pde.elements.lev(idx,:)-1; % NOTE : remove the 1 per note below
            thisElemPosVec = pde.elements.pos(idx,:)-1; % NOTE : remove the 1 per note below
            
            for d=1:nDims
                
                %%
                % Now add two daughter nodes in each dimension
                
                for d2=1:nDims
                    
                    %%
                    % TODO : This refinement is not sparse. Need to add the
                    % sparse grid selection rule here also.
                    
                    %%
                    % First daughter
                    
                    newElemLevVec = thisElemLevVec;
                    newElemPosVec = thisElemPosVec;
                    
                    newElemLevVec(d2) = thisElemLevVec(d)+1;
                    newElemPosVec(d2) = thisElemPosVec(d)*2-1;
                    
                    element_idx = lev_cell_to_element_index(pde,newElemLevVec,newElemPosVec);
                    pde.elementsIDX(nElements+cnt) = element_idx; % Extend element list
                    i1 = (nElements+cnt-1)*elementDOF-1; % Get the start and end global row indices of the new element
                    i2 = (nElements+cnt)*elementDOF;
                    fval(i1:i2) = newElementVal; % Extend coefficient list with near zero magnitude (ideally would be zero)
                    
                    %                     for d=1:nDims
                    %                         pde.elements.coords{d}.lev (element_idx) = newElemLevVec(d)+1; % NOTE : have to start lev  index from 1 for sparse storage
                    %                         pde.elements.coords{d}.cell(element_idx) = newElemPosVec(d)+1; % NOTE : have to start cell index from 1 for sparse storage
                    %                     end
                    %                     for d=1:nDims
                    pde.elements.lev(element_idx,:) = newElemLevVec+1; % NOTE : have to start lev  index from 1 for sparse storage
                    pde.elements.pos(element_idx,:) = newElemPosVec+1; % NOTE : have to start cell index from 1 for sparse storage
                    %                     end
                    
                    cnt = cnt + 1;
                    
                    %%
                    % Second daughter
                    
                    newElemLevVec = thisElemLevVec;
                    newElemPosVec = thisElemPosVec;
                    
                    newElemLevVec(d2) = thisElemLevVec(d)+1;
                    newElemPosVec(d2) = thisElemPosVec(d)*2;
                    
                    element_idx = lev_cell_to_element_index(pde,newElemLevVec,newElemPosVec);
                    pde.elementsIDX(nElements+cnt) = element_idx; % Extend element list
                    i1 = (nElements+cnt-1)*elementDOF-1; % Get the start and end global row indices of the new element
                    i2 = (nElements+cnt)*elementDOF;
                    fval(i1:i2) = newElementVal; % Extend coefficient list with near zero magnitude (ideally would be zero)
                    
                    %                     for d=1:nDims
                    %                         pde.elements.coords{d}.lev (element_idx) = newElemLevVec(d)+1; % NOTE : have to start lev  index from 1 for sparse storage
                    %                         pde.elements.coords{d}.cell(element_idx) = newElemPosVec(d)+1; % NOTE : have to start cell index from 1 for sparse storage
                    %                     end
                    %                     for d=1:nDims
                    %                         pde.elements.lev(element_idx) = newElemLevVec+1; % NOTE : have to start lev  index from 1 for sparse storage
                    %                         pde.elements.pos(element_idx) = newElemPosVec+1; % NOTE : have to start cell index from 1 for sparse storage
                    %                     end
                    
                    pde.elements.lev(element_idx,:) = newElemLevVec+1; % NOTE : have to start lev  index from 1 for sparse storage
                    pde.elements.pos(element_idx,:) = newElemPosVec+1; % NOTE : have to start cell index from 1 for sparse storage
                    
                    cnt = cnt + 1;
                    
                end
                
            end
            
        else
            
            fprintf('leaf node but no refinement, fval=%f\n',fval(gidx))
        end
        
    else
        
        disp('internal node');
        
    end
    
end

assert(numel(fval)==numel(pde.elementsIDX)*elementDOF);

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

end