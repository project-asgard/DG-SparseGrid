function [hash_table,A_data,f_new] = addHierLimitElements(pde,opts,hash_table,A_data,f,M,m)
%Q represents the polynomial coefficents in the realspace tensor DG space.

persistent Ix FG2DG 

assert(numel(pde.dimensions) == 2);
assert(opts.deg == 2); %Stick with linear polynomials for now

pde_lev_vec = [pde.dimensions{1}.lev,pde.dimensions{2}.lev];
pde_lev_max = max(pde_lev_vec);
deg = opts.deg;

pos_tol = 1e-14;

%% Create kron 2 DG transformation matrix.  This is independent of adaptivity
if isempty(FG2DG)
    FG2DG = cell(max(pde_lev_vec));
    Ix = cell(max(pde_lev_vec),1);
    
    
    for l=1:pde_lev_vec(1) %Level 0 is always in the space.
        for k=1:pde_lev_vec(2)
            lev_vec = [l k];
            
            dx = (pde.dimensions{1}.max-pde.dimensions{1}.min)/2^lev_vec(1);
            dv = (pde.dimensions{2}.max-pde.dimensions{2}.min)/2^lev_vec(2);
        
            %Build restriction matrices to lower level space
            R_x = speye(2^lev_vec(1)*deg,2^pde_lev_vec(1)*deg);
            R_v = speye(2^lev_vec(2)*deg,2^pde_lev_vec(2)*deg);
            R   = kron(R_x,R_v);
            FMWT_x = OperatorTwoScale_wavelet2(deg,lev_vec(1));
            FMWT_v = OperatorTwoScale_wavelet2(deg,lev_vec(2));
            FMWT_2D = kron(FMWT_x',FMWT_v');
            fg2dg    = kronrealspace2DtoDG(lev_vec,deg,deg)*FMWT_2D;
            FG2DG{l,k} = fg2dg*R;
            


        end
        %Lastly construst wavelet passthrough for representing local
        %elements in the hierarchical wavelet space
        [Ix{l},~] = find(OperatorTwoScale_wavelet2(1,l));
    end 
    B0 = kronrealspace2DtoDG(pde_lev_vec,deg,1)*FMWT_2D;
end



%% Map SG to FG space
[perm,~,pvec] = sg_to_fg_mapping_2d(pde,opts,A_data);

f_FG = zeros(size(pvec));
f_FG(pvec) = f(perm(pvec));


%% Create hierarchical order 

[A,B] = meshgrid(1:pde_lev_vec(1),1:pde_lev_vec(2));
C = cat(2,A',B');
C = reshape(C,[],2);

%Sort method 1
% [s_sort,I] = sort(sum(C,2));
% pos_lev_pde = C(I,:);
% pos_lev_pde(s_sort <= max(pde_lev_vec),:) = [];

%Sort method 2
[~,I] = sort(sum(C,2)*(pde_lev_max+1)+abs(diff(C,1,2)));
pos_lev_pde = C(I,:);
pos_lev_pde(max(pos_lev_pde,[],2) ~= pde_lev_max,:) = [];
pos_lev_pde(min(pos_lev_pde,[],2) < 4,:) = [];

x_F = pde.dimensions{1}.min:(pde.dimensions{1}.max-pde.dimensions{1}.min)/2^pde_lev_vec(1):pde.dimensions{1}.max;
v_F = pde.dimensions{2}.min:(pde.dimensions{2}.max-pde.dimensions{2}.min)/2^pde_lev_vec(2):pde.dimensions{2}.max;


%% Loop over hierarchy
for l=1:size(pos_lev_pde,1)
    %Get current level
    lev_vec = pos_lev_pde(l,:);
    fprintf('Lev = [%d,%d].\n',lev_vec(1),lev_vec(2));
    
    %Transform to DG realspace
    Q = FG2DG{lev_vec(1),lev_vec(2)}*f_FG;
    
    %Need higher-level perturbation that is not affected by the limiter
    perturb = f_FG - FG2DG{lev_vec(1),lev_vec(2)}'*Q;
    
    %Useful for debugging 
    x = [pde.dimensions{1}.min:(pde.dimensions{1}.max-pde.dimensions{1}.min)/2^lev_vec(1):pde.dimensions{1}.max];
    v = [pde.dimensions{2}.min:(pde.dimensions{2}.max-pde.dimensions{2}.min)/2^lev_vec(2):pde.dimensions{2}.max];
    
    figure(19);
    subplot(2,2,1); plotVec(x_F,v_F,1,FG2DG{pde_lev_vec(1),pde_lev_vec(2)}*f_FG)
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Fine Grid Function';title({mytitle,old_title});
    subplot(2,2,2); plotVec(x_F,v_F,0,B0*f_FG)
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Fine Grid Cell Average';title({mytitle,old_title});
    subplot(2,2,3); plotVec(x,v,1,Q)
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Projected Function';title({mytitle,old_title});
    subplot(2,2,4); plotVec(x,v,0,Q(1:4:end))
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Projected Cell Average';title({mytitle,old_title});
    sgtitle(sprintf('Before Limiting on lev = [%d,%d]',lev_vec(1),lev_vec(2)));
    %Add elements if negative on this level
    [limited,hash_table,Q] = hierLimit(pde,lev_vec,pde_lev_vec,opts,hash_table,Q,Ix{lev_vec(1)},Ix{lev_vec(2)},M,m,pos_tol);
    if limited
        %Update f
        f_FG = FG2DG{lev_vec(1),lev_vec(2)}'*Q + perturb;
    end
    figure(20);
    subplot(2,2,1); plotVec(x_F,v_F,1,FG2DG{pde_lev_vec(1),pde_lev_vec(2)}*f_FG)
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Fine Grid Function';title({mytitle,old_title});
    subplot(2,2,2); plotVec(x_F,v_F,0,B0*f_FG)
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Fine Grid Cell Average';title({mytitle,old_title});
    subplot(2,2,3); plotVec(x,v,1,Q)
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Projected Function';title({mytitle,old_title});
    subplot(2,2,4); plotVec(x,v,0,Q(1:4:end))
    h1 = get(gca,'title');old_title = get(h1,'string');mytitle = 'Projected Cell Average';title({mytitle,old_title});
    sgtitle(sprintf('After Limiting on lev = [%d,%d]',lev_vec(1),lev_vec(2)));
end





%% Modify coordinates to produce limited f

%Get A_data for new hash table
A_data = global_matrix(pde,opts,hash_table);

%Get SG to FG transformation
[~,iperm,~] = sg_to_fg_mapping_2d(pde,opts,A_data);

%Now to adpative space
f_new = f_FG(iperm);
end

function [limited,hash_new,lQ] = hierLimit(pde,lev_vec,pde_lev_vec,opts,hash_table,Q,Ix,Iv,M,m,tol)

%%
limiter = 'Pos';


n_v = uint64(2^lev_vec(2));
lev_x = lev_vec(1);
lev_v = lev_vec(2);
max_lev = max(pde_lev_vec);
num_dims = numel(lev_vec);

%%% Limiter data  %%%%%%%%%%%%%%%%%%
dx = (pde.dimensions{1}.max-pde.dimensions{1}.min)/2^lev_x;
dv = (pde.dimensions{2}.max-pde.dimensions{2}.min)/2^lev_v;
%Transformation from DG basis to Limiter basis
L = diag([sqrt(1/(dx*dv)),sqrt(3/(dx*dv)),sqrt(3/(dx*dv)),3*sqrt(1/(dx*dv))]);
%Transformation from Limiter basis to DG basis
Linv = inv(L);
%Evaluate the function at the verticies of the rectangle
evalmat = [1 1 1 1;1 1 -1 -1;1 -1 1 -1; 1 -1 -1 1];
%Vector for cell average
e1 = zeros(4,1); e1(1) = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

I = [];
lQ = Q; %%Limited Q
limited = true;
%% Limit solution %%%%%%%%%%%%%%%%%%
num_x = 2^lev_x; num_v = 2^lev_v; k = 1;

if strcmp(limiter,'TVD') %%%%% TVD Limiter %%%%%%%%%%%%%%%%
    for i=1:num_x
        for j=1:num_v
            %Convert to proper basis
            %Get local index
            start = ((i-1)*num_v+(j-1))*(k+1)^2;
            curr = L*Q(start+1:start+(k+1)^2);
            curr(1) = max([curr(1),0]);

            up = zeros(4,1);
            dn = zeros(4,1);
            lt = zeros(4,1);
            rt = zeros(4,1);

            %Get elements to up, down, left, and right of u
            if j < num_v
                idx = start + (k+1)^2;
                up = L*Q(idx+1:idx+(k+1)^2);
                up(1) = max([up(1),0]);
            end

            if j > 1
                idx = start - (k+1)^2;
                dn = L*Q(idx+1:idx+(k+1)^2);
                dn(1) = max([dn(1),0]);
            end

            if i > 1
                idx = start - num_v*(k+1)^2;
                lt = L*Q(idx+1:idx+(k+1)^2);
                lt(1) = max([lt(1),0]);
            end

            if i < num_x
                idx = start + num_v*(k+1)^2;
                rt = L*Q(idx+1:idx+(k+1)^2);
                rt(1) = max([rt(1),0]);
            end

            %if abs(curr(3)) >= M*dx^2 || abs(curr(2)) >= M*dv^2

                uy = minmod(curr(2),(up(1)-curr(1))/2,(curr(1)-dn(1))/2);
                ux = minmod(curr(3),(rt(1)-curr(1))/2,(curr(1)-lt(1))/2);
            if norm([uy,ux]-curr(2:3)) > 1e-15
                %fprintf('LIMIT:    Limiting on (%3d,%3d)\n',i,j);
                %curr(2) = 0;
                %curr(3) = 0;
                curr(2) = uy;
                curr(3) = ux;
                curr(4) = 0;

                %if ~quiet
                %    fprintf('Limit in x for cell with center (%f,%f)\n',0.5*(x(i)+x(i+1)),0.5*(v(j)+v(j+1)));
                %end
                %if k == 2
                %    ru(start+5:start+9) = zeros(5,1);
                %end

                %Add this element to the list of elements to add if needbe
                I = [I;uint64((i-1)*num_v+(j-1)+1)];
            end

            lQ(start+1:start+(k+1)^2) = Linv*curr;
        end
    end
elseif strcmp(limiter,'Pos') %%%%% Maximum Principle Limiter %%%%%%%%%%%%%
    for i=1:num_x
        for j=1:num_v
            start = ((i-1)*num_v+(j-1))*(k+1)^2;
            curr = L*Q(start+1:start+(k+1)^2);

            %Need max and min for must happen at verticies of rectangle
            uval = evalmat*curr;
            u_M = max(uval);
            u_m = min(uval);
            %Limiter modifier
            if u_M > M+tol || u_m < m-tol
                %curr(1) = max([curr(1),0]);
                theta = min(abs([1,((M+tol)-curr(1))/(u_M-curr(1)),((m-tol)-curr(1))/(u_m-curr(1))]));
                %if abs(theta-1) > 5e-15 %Need to limit
                %fprintf('Limiting on cell (%d,%d) with err %e\n',i,j,abs(theta-1));
                %fprintf('--> Theta_vec = [%f,%f,%f]\n',abs([1,(M-curr(1))/(u_M-curr(1)),(m-curr(1))/(u_m-curr(1))]));
                I = [I;uint64((i-1)*num_v+(j-1)+1)];
                u_avg = curr(1)*e1;
                curr = theta*(curr-u_avg)+u_avg;
                lQ(start+1:start+(k+1)^2) = Linv*curr;
            end
            
        end
    end
end

if isempty(I)
    limited = false;
    lQ = Q;
    return
end
    

%% Get elements to represent 1 on that cell
%Get x and v indicies of I
x_idx = idivide(I-1,n_v) + 1;
v_idx = mod(I-1,n_v) + 1;

pre_elements_to_add = zeros(numel(x_idx)*(lev_x+1)*(lev_v+1),2);
%Need to go from realspace to hierarchical space
for i=1:numel(x_idx)
    x_hier = Ix((x_idx(i)-1)*(lev_x+1)+1:x_idx(i)*(lev_x+1));
    v_hier = Iv((v_idx(i)-1)*(lev_v+1)+1:v_idx(i)*(lev_v+1));
    
    %Matlab magic to do list all combinations of these two vectors from
    % -------
    % https://www.mathworks.com/matlabcentral/answers/98191-how-can-i-obtain-all-possible-combinations-of-given-vectors-in-matlab
    % -------
    [A,B] = meshgrid(x_hier,v_hier);
    C = cat(2,A',B');
    pre_elements_to_add( (i-1)*(lev_x+1)*(lev_v+1)+1:i*(lev_x+1)*(lev_v+1),:) = reshape(C,[],2);
end

%Get global md index of each element
idx_vec = pre_elements_to_add(:,1) + (pre_elements_to_add(:,2)-1)*2^max_lev;

%% Add elements to a new hash table

%Get lev and position of each element
hash_new = hash_table;

num_elements = numel(hash_table.elements_idx);
num_elements_added = 0;
for i=1:numel(idx_vec)
    idx = idx_vec(i);
    
    if hash_new.elements.type(idx) == 0 % element not already enabled
            
        num_elements_added = num_elements_added + 1;
        position_in_elements_idx = num_elements+num_elements_added;
        hash_new.elements_idx(position_in_elements_idx) = idx; % Extend element list
        
        [lev_vec, pos_vec] = md_idx_to_lev_pos(num_dims, opts.max_lev, idx);

        hash_new.elements.lev_p1(idx,:) = lev_vec+1; % NOTE : have to start lev  index from 1 for sparse storage
        hash_new.elements.pos_p1(idx,:) = pos_vec+1; % NOTE : have to start cell index from 1 for sparse storage
        hash_new.elements.type(idx) = 1;

    end
    
end


end




function z = minmod(a,b,c)
    if abs(sign(a)+sign(b)+sign(c)) == 3 %Equiv to sign(a) = sign(b) = sign(c)
        s = sign(a);
        z = s*min(abs([a,b,c]));
    else
        z = 0;
    end
end

