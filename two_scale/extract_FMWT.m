function FMWT = extract_FMWT( pde, dimIdx, FMWT )
% FMWT = extract_FMWT( pde, dimIdx, FMWT )
% 
% --------------------------------------------------------------------------
% Now extract out only the elements we need (i.e., for the case on non-full
% levels during adativity)
% --------------------------------------------------------------------------
do_output = 1;

    N = numel(pde.elementsIDX);
    NF = numel(FMWT(:,1));
    keep = zeros(NF,1);
    
%     thisDimLev = pde.elements.lev;
%     thisDimPos = pde.elements.pos;
    
    for n=1:N
        
        idx = pde.elementsIDX(n);
        
%         thisLev = full(thisDimLev(idx)-1);
%         thisPos = full(thisDimPos(idx)-1);
        thisLev = full(pde.elements.lev(idx,dimIdx)-1);
        thisPos = full(pde.elements.pos(idx,dimIdx)-1);
        
        idx1D = lev_cell_to_singleD_index(thisLev,thisPos);
        if (do_output),
          fprintf('lev: %i, pos: %i, idx: %i\n',thisLev,thisPos,idx1D);
        end;
        keep((idx1D-1)*deg+1:idx1D*deg) = 1;
    end
    
    jdx = find(keep);
    FMWT = FMWT(jdx,jdx);
end
