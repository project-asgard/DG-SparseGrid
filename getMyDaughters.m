function [newElemLevVecs,newElemPosVecs,cnt] = getMyDaughters(pde,idx)

nDims = numel(pde.dimensions);

debug = 1;

%%
% Get this coordinate vector
% levVec = [lev1,lev2,...,levD]
% posVec = [pos1,pos2,...,posD]

thisElemLevVec = pde.elements.lev_p1(idx,:)-1; % NOTE : remove the 1 per note below
thisElemPosVec = pde.elements.pos_p1(idx,:)-1; % NOTE : remove the 1 per note below

if thisElemLevVec(1) == 2
    disp('should never get here');
end

newElemLevVecs = [];
newElemPosVecs = [];

%%
% Generate list of new elements which satisfy selection rule

cnt = 0;
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
        
        if sum(newElemLevVec)<=thisElemLevVec(d)+1 & newElemLevVec(d2)<=pde.maxLev % Sparse grid selection rule AND max depth check
            
            newElemLevVecs(cnt+1,:) = newElemLevVec;
            newElemPosVecs(cnt+1,:) = newElemPosVec; % Assumes pos starts at 0
            
            assert(newElemPosVecs(cnt+1,d2) >= 0);
            assert(newElemLevVecs(cnt+1,d2) >= 0);
            
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
        
        if sum(newElemLevVec)<=thisElemLevVec(d)+1 & newElemLevVec(d2)<=pde.maxLev % Sparse grid selection rule AND max depth check
            
            newElemLevVecs(cnt+1,:) = newElemLevVec;
            newElemPosVecs(cnt+1,:) = newElemPosVec; % Assumes pos starts at 0
            
            assert(newElemPosVecs(cnt+1,d2) >= 0);
            assert(newElemLevVecs(cnt+1,d2) >= 0);
            
            cnt = cnt + 1;
            
            pde.elements.node_type(idx) = 1; % Now that this element has been refined it is no longer a leaf.
            
        else
            
            disp('element not added because it did not obey sparse selection rule');
            
            
        end
        
    end
    
end

end