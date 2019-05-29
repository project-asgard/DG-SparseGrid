function fval_r = real_space_solution_at_sparse_grid(pde,wavelet_coeffs)

nElem = numel(pde.elementsIDX);
nDims = numel(pde.dimensions);
nDOF = numel(wavelet_coeffs);
deg = pde.dimensions{1}.deg;

assert(nDOF == nElem*deg^nDims);

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

allCoords(:,1) = [-1:0.005:+1];

nPts = numel(allCoords(:,1));

fval_r = zeros(nPts,1);

%%
% For the desired set of grid points, constuct (for each dimension) the M
% transform ...
%
% M * f_legendre => f_realspace


MM = sparse(nPts,nDOF);
for d=1:nDims
    
    dim         = pde.dimensions{d};
    max_lev      = dim.lev;
    nLegendre   = 2^max_lev;
    
    %%
    % Get the scale_fac for this level
    
    dMin = dim.domainMin;
    dMax = dim.domainMax;
    h = (dMax-dMin)/nLegendre;
    scale_fac = sqrt(1/h);
    
    x = allCoords(:,d);
    
    clear MM;
    for l=1:nLegendre
        
        %%
        % Get the normalized to [-1,+1] coordinates for this element
        
        xMin = dMin + (l-1)*h;
        xMax = xMin + h;
        x_norm = (x - xMin)/(xMax-xMin)*2-1;
        
        %%
        % Evaluate the legendre functions at these coordinates with
        % P(x)<-1 == 0 & P(x)>+1 == 0
        
        p_val = lin_legendre(x_norm,deg) * scale_fac;
        
        %%
        % Get the lev, pos and range of this element
        
        %         idx = pde.elementsIDX(l);
        %         [elem_coord_D,elem_coord_L_D,elem_coord_R_D] = getMyRealSpaceCoord(pde,idx);
        %         lev_D = pde.elements.lev_p1(idx,:)-1;
        %         lev = lev_D(d);
        
        %%
        % Plot realspace basis functions
        
        figure(88)
%         cla
        for dd=1:deg
            plot(x,1./(1./(p_val(:,dd))),'-o');
            hold on
            ylim([-10,10]);
        end
%         hold off
        
        j1 = deg*(l-1)+1;
        j2 = deg*(l);
        
        MM(:,j1:j2) = p_val;
        
    end
    
    M_D{d} = MM;
    
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
        
%         IDlev = pde.elements.lev_p1(pde.elementsIDX(elem),d)-1;
%         IDpos = pde.elements.pos_p1(pde.elementsIDX(elem),d)-1;
%         IDe = lev_cell_to_singleD_index(IDlev,IDpos);
%         
%         assert(IDe==elem);
        
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

figure()
x = allCoords(:,1);
[X,I] = sort(x);
plot(X,fval_r(I),'o')

end