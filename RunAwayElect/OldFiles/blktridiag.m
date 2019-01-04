function A = blktridiag(Amd,Asub,Asup,n)

% get block sizes, check for consistency
[p,q] = size(Amd);


% scalar inputs?
% since p and q are integers...
if (p*q)==1
    if n==1
        A = Amd;
    else
        % faster as Jos points out
        A = spdiags(repmat([Asub Amd Asup],n,1),-1:1,n,n);
    end
    % no need to go any farther
    return
end

% use sparse. the main diagonal elements of each array are...
v = repmat(Amd(:),n,1);
% then the sub and super diagonal blocks.
if n>1
    % sub-diagonal
    v=[v;repmat(Asub(:),n-1,1)];
    
    % super-diagonal
    v=[v;repmat(Asup(:),n-1,1)];
end


% now generate the index arrays. first the main diagonal
[ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-1);
rind = 1+ind1(:)+p*ind3(:);
cind = 1+ind2(:)+q*ind3(:);
% then the sub and super diagonal blocks.
if n>1
    % sub-diagonal
    [ind1,ind2,ind3]=ndgrid(0:p-1,0:q-1,0:n-2);
    rind = [rind;1+p+ind1(:)+p*ind3(:)];
    cind = [cind;1+ind2(:)+q*ind3(:)];
    
    % super-diagonal
    rind = [rind;1+ind1(:)+p*ind3(:)];
    cind = [cind;1+q+ind2(:)+q*ind3(:)];
end

% build the final array all in one call to sparse
A = sparse(rind,cind,v,n*p,n*q);

% ====================================
% end mainline
% ====================================