%%
% Now extract out only the elements we need (i.e., for the case on non-full
% levels during adativity)

extract = 0;
if extract
    N = numel(pde.elementsIDX);
    NF = numel(FMWT(:,1));
    keep = zeros(NF,1);
    
    thisDimLev = pde.elements.coords{dimIdx}.lev;
    thisDimPos = pde.elements.coords{dimIdx}.cell;
    
    for n=1:N
        
        idx = pde.elementsIDX(n);
        
        thisLev = full(thisDimLev(idx)-1);
        thisPos = full(thisDimPos(idx)-1);
        
        idx1D = lev_cell_to_singleD_index(thisLev,thisPos);
        fprintf('lev: %i, pos: %i, idx: %i\n',thisLev,thisPos,idx1D);
        keep((idx1D-1)*deg+1:idx1D*deg) = 1;
    end
    
    FMWT = FMWT(find(keep),find(keep));
end