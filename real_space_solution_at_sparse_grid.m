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
% For the desired set of grid points, constuct (for each dimension) the M
% transform ...
%
% M * f_legendre => f_realspace

for d=1:nDims
    
    dim = pde.dimensions{d};
    
    clear M;
    for elem=1:nElem
        
        %%
        % Get the lev, pos and range of this element
        
        idx = pde.elementsIDX(elem);
        [elem_coord_D,elem_coord_L_D,elem_coord_R_D] = getMyRealSpaceCoord(pde,idx);
        lev_D = pde.elements.lev_p1(idx,:)-1;
        lev = lev_D(d);
        
        %%
        % Get the scale_fac for this element
        
        xMin = elem_coord_L_D(d); % left side of element
        xMax = elem_coord_R_D(d); % right side of element
        h = (xMax-xMin)/(2^lev);
        scale_fac = sqrt(1/h);
        
        %%
        % Get the normalized to [-1,+1] coordinates for this element
        
        x = allCoords(:,d);
        x_norm = (x - xMin)/(xMax-xMin)*2-1;
        
        %%
        % Evaluate the legendre functions at these coordinates with
        % P(x)<-1 == 0 & P(x)>+1 == 0
        
        p_val = lin_legendre(x_norm,deg) * scale_fac;
        
        %         %%
        %         % Plot realspace basis functions
        %
        %         figure(88)
        %         cla
        %         for dd=1:deg
        %             plot(x_norm,p_val(:,dd));
        %             hold on
        %         end
        %         hold off
        
        i1 = deg*(elem-1)+1;
        i2 = deg*(elem);
        
        M(:,i1:i2) = p_val;
        
    end
    
    M_D{d} = M;
    
end

%%
% Evaluate M * F' * f_wavelet = f_realspace
%
% Other notes ...
%
% M  * f_legendre => f_realspace
% F' * f_wavelet  => f_legendre

for elem=1:nElem
    
    clear kronMatList;
    for d=1:nDims
        
        i1 = deg*(elem-1)+1;
        i2 = deg*(elem);
        
        FMWT_T = dim.FMWT';
        F_T = FMWT_T(i1:i2,i1:i2); % note the transpose
        
        M = M_D{d}(:,i1:i2);
        
        kron_mat_list{d} = M * F_T; % Build kron list
        
    end
    
    fval_r =  fval_r + kron_multd(nDims,kron_mat_list,wavelet_coeffs);
    
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