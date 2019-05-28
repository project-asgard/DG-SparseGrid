function fval_r = real_space_solution_at_sparse_grid(pde,wavelet_coeffs)

nElem = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);
deg = pde.dimensions{1}.deg;

%%
% Get center coordinates for all elements

% allCoords = zeros(N,nDims);
% for n=1:N
%     idx = pde.elementsIDX(n);
%     allCoords(n,:) = getMyRealSpaceCoord(pde,idx);
% end

%%
% For now just use a simple, ordered set of coordinates (the above would
% give the centers of all the elements)

allCoords(:,1) = [-1:0.01:+1];

nPts = numel(allCoords(:,1));

fval_r = zeros(nPts,1);

%%
% Loop over elements, and evaluate which other elements each contributes to
% and evaluate and add its contribution.

for elem=1:nElem
    
    %%
    % Get this element coordinate vector
    
    idx = pde.elementsIDX(elem);
    
    [myCoord,myCoordL,myCoordR] = getMyRealSpaceCoord(pde,idx);
    lev_D = pde.elements.lev_p1(idx,:)-1;
    pos_D = pde.elements.pos_p1(idx,:)-1;
    
    clear kronMatList;
    for d=1:nDims
             
        dim = pde.dimensions{d};
        lev = lev_D(d);
        pos = pos_D(d);
        
        %%
        % Get this elements piece of 1D FMWT (legendre -> wavelet)
        
        idx_1D = lev_cell_to_singleD_index(lev,pos);
        i1 = (idx_1D-1)*deg+1;
        i2 = idx_1D*deg;
        
        FMWT = dim.FMWT;
        F_T = FMWT(i1:i2,i1:i2)'; % note the transpose
        
        %%
        % Get this elements scale factor
        
        dMin = dim.domainMin;
        dMax = dim.domainMax;
        h = (dMax-dMin)/(2^lev);
        scale_fac = sqrt(1/h);
        
        %%
        % Map desired coordinates to this elements 1D [-1,+1] range
        
        xMin = myCoordL(d); % left side of element
        xMax = myCoordR(d); % right side of element
        
        x = allCoords(:,d);
        x_norm = (x - xMin)/(xMax-xMin)*2-1;
        
        %%
        % Get legenre functions at these locations (legendre -> realspace)
        
        M = lin_legendre(x_norm,deg) * scale_fac;
        
        %%
        % F_T = wavelet  -> legendre
        % M   = legendre -> real
        
        kron_mat_list{d} = M * F_T; % Build kron list
        
    end
    
    fval_r = fval_r + kron_multd(nDims,kron_mat_list,wavelet_coeffs);
    
    %     figure(99);
    %     hold on
    %     x = allCoords(:,1);
    %     [X,I] = sort(x);
    %     plot(X,fval_r(I),'-o')
    
end

%%
% Test plot

x = allCoords(:,1);
[X,I] = sort(x);
plot(X,fval_r(I),'o')

end