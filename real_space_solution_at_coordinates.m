function fval_r = real_space_solution_at_coordinates(pde,wavelet_coeffs,coordinates)

num_elem = numel(pde.elementsIDX);
num_dims = numel(pde.dimensions);
num_DOF  = numel(wavelet_coeffs);
num_pts  = numel(coordinates(:,1));

deg      = pde.dimensions{1}.deg;

assert(num_DOF == num_elem * deg^num_dims);

%%
% For the desired set of grid points, constuct (for each dimension) the M
% transform ...
%
% M * f_legendre -> f_realspace

MM = sparse(num_pts,num_DOF);
for d=1:num_dims
     
    dim             = pde.dimensions{d};
    max_lev         = dim.lev;
    num_legendre    = 2^max_lev;
     
    %%
    % Get the scale_fac for this level
    
    d_min     = dim.min;
    d_max     = dim.max;
    h         = (d_max-d_min) / num_legendre;
    scale_fac = sqrt(1/h);
    
    x = coordinates(:,d);
    
    clear MM;
    for l=1:num_legendre
        
        %%
        % Get the normalized to [-1,+1] coordinates for this element
        
        x_min = d_min + (l-1)*h;
        x_max = x_min + h;
        x_norm = (x - x_min)/(x_max-x_min)*2-1;
        
%         %%
%         % Check for evaluations at multi-valued locations
%         
%         iiMulti = find(x_norm == -1 | x_norm == +1);
%         TODO
        
        
        %%
        % Evaluate the legendre functions at these coordinates with
        % P(x)<-1 == 0 & P(x)>+1 == 0
        
        p_val = lin_legendre(x_norm,deg) * scale_fac;
        
        j1 = deg*(l-1)+1;
        j2 = deg*(l);
        
        MM(:,j1:j2) = p_val;
        
    end
    
    M_D{d} = MM;
    
end


%%
% Evaluate the transform from wavelet-> realspace, i.e., 
%
% since
%
% F = FMWT
%
% F' * f_wavelet  -> f_legendre
%
% M  * f_legendre -> f_realspace
%
% so
%
% M * F' * f_wavelet -> f_realspace

for d=1:num_dims
    right_trans = 'RT';
    %MF_T{d} = MM * dim.FMWT';
    MF_T{d} = apply_FMWT_blocks(pde.dimensions{d}.lev, FMWT_blocks, MM, right_trans);
end

%%
% Perform kron product to combine dimensions

fval_r = zeros(num_pts,1);

for elem=1:num_elem
    
    clear kronMatList;
    for d=1:num_dims
             
        %%
        % the 1D index requires a mapping because of adaptivity
        
        lev = pde.elements.lev_p1(pde.elementsIDX(elem),d)-1;
        pos = pde.elements.pos_p1(pde.elementsIDX(elem),d)-1;
        idx_1D = lev_cell_to_singleD_index(lev,pos);
        
%         assert(idx_1D==elem);
        
        i1 = deg*(idx_1D-1)+1;
        i2 = deg*(idx_1D);
        
        kron_mat_list{d} = MF_T{d}(:,i1:i2);
        
    end
    
    element_ii = deg^num_dims*(elem-1)+1:deg^num_dims*elem;

    fval_r = fval_r + kron_multd(num_dims, kron_mat_list, wavelet_coeffs(element_ii) );
    
end

%%
% Test plot

do_plot = 0;
if do_plot
    figure(2000)
    if num_dims == 1
        x = coordinates(:,1);
        plot(x,fval_r,'o');
    elseif num_dims == 2
    end
end

end
