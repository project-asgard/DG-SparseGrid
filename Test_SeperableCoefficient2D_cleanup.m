% ------------------------------------------------
% Implementation of Stationary Advection equation
% in 2D:
%
%  div (A u) = f in domain Omega
%         u  = g on the in flow boundary
% *********************
% Here A is a vector A = <k,h>^T 
% and k(x,y) and h(x,y) are non-seperable coefficients 
% ------------------------------------------------

clear
close all

% ==============================
% Parameters for using this code
% ==============================
use_kronmult = 1;


lev = 5;
deg = 3;
xMin = 0;xMax = 1;
yMin = 0;yMax = 1; % For simplicity, I assume Omega_y = Omega_x

% Test 
alpha = 1;
beta = 1;
kcoef = @(x,y)sin(beta*x+alpha*y); 
hcoef = @(x,y)sin(beta*x+alpha*y); 
exactu = @(x,y)sin(pi*x).*sin(pi*y);
FunF = @(x,y)sin(pi*x).*sin(pi*y).*(beta + alpha).*cos(alpha*y + beta*x) + sin(alpha*y + beta*x)*pi.*(cos(pi*x).*sin(pi*y) + cos(pi*y).*sin(pi*x));




% =============================================
% Generate the kcoef and hcoef corresponding to 
% A = <k,h>^T 
% =============================================

% Setup jacobi of variable x and define coeff_mat
N = 2^(lev);
h = (xMax-xMin) / N;
dof_1D = deg * N;
Jacobi = h/2;

%%
%  Get quadrature points and weights.
%  quad_x(:) is quad_num by 1
%  quad_w(:) is quad_num by 1
[quad_x,quad_w] = lgwt(default_quad_number(deg),-1,1);

FMWT_COMP = OperatorTwoScale(deg,lev);

%  Compute the trace values (values at the left and right of each element for all k)
%  p_L(:) is 1 by deg
%  p_R(:) is 1 by deg
% The trace values
p_L = lin_legendre(-1,deg) * 1/sqrt(h);
p_R = lin_legendre(+1,deg) * 1/sqrt(h);

%%
%  Get the basis functions and derivatives for all k
%  p_val(:,:) is quad_num by deg
%  Dp_val(:,:) is quad_num by deg
p_val  = lin_legendre(quad_x,deg)  * 1/sqrt(h);
Dp_val = lin_dlegendre(quad_x,deg) * 1/sqrt(h) * 2/h;

% =======================================================
% Obtain the k_vect and h_vect
% k = sum k_vect*phi(x)*psi(y)
% h = sum h_vect*phi(x)*psi(y)
% Note: I am using full grid + non-multiwavelet for these
% two coefficients
% =======================================================
k_vect = zeros((N*deg)^2,1);
h_vect = zeros((N*deg)^2,1);
for i = 0:N-1
    
    for j=0:N-1
        
        yL = xMin + j*h;
        yR = yL + h;
        
        quad_xi_y = (((quad_x+1)/2+j)*h+xMin);
        
        num_cell = N*(j)+i+1;
        
        xL = xMin + i*h;
        xR = xL + h;
        
        % Map quadrature points from [-1,1] to physical domain of this i element
        
        quad_xi_x = (((quad_x+1)/2+i)*h+xMin);
        
        [xx,yy] = meshgrid(quad_xi_x,quad_xi_y);
        tmp_k = kcoef(xx,yy);
        tmp_h = hcoef(xx,yy);
        
        p_val2D = kron(p_val,p_val);
        quad_w_2D = kron(quad_w,quad_w);
        
        val_mass_k = p_val2D' * (tmp_k(:) .* quad_w_2D) * Jacobi*Jacobi;
        val_mass_h = p_val2D' * (tmp_h(:) .* quad_w_2D) * Jacobi*Jacobi;
        
        for k = 1:deg           
              Ind(deg*(k-1)+[1:deg]) = i*(N*deg^2)+j*deg+deg*N*(k-1)+[1:deg];
        end
        
        % for projection of k-coefficient function
        k_vect(Ind) = val_mass_k;
        h_vect(Ind) = val_mass_h;
        
    end
end



% ===========================================
% Generate matrices for Mass and Grad
% Mass: int phi_i*phi_j*phi_k dx
% Grad: 
% (grad k*u,v) = -(k*u,grad v) + <hat{k*u},v>
% Here k will be replaced by the DG basis
% 
% ===========================================
tmp1D = zeros(N*deg,N*deg,N*deg);

BCL = 'D';
BCR = 'N';

FluxVal = 1; % This is central flux
FCL = 1; FCR = 1;
for i = 0:N-1
    c = deg*i + [1:deg];
    for l = 1:deg
        % (dv/dx*k,u)
        val_grad  = -Dp_val'* (p_val .* p_val(:,l) .* quad_w) * Jacobi;
        tmp1D(c,c,deg*i+l) = tmp1D(c,c,deg*i+l) + val_grad;
        %%
        % - <funcCoef*{q},p>
        %----------------------------------------------
        % Numerical Flux is defined as
        % Flux = {{f}} + C/2*[[u]]
        %      = ( f_L + f_R )/2 + FunCoef*( u_R - u_L )/2
        % [[v]] = v_R - v_L
         
        FCL = 1;
        FCR = 1;

        TraVal = [...
            (-p_L') * FCL * p_R/2.*p_L(:,l) + FluxVal * abs(FCL)/2 * (-p_L') *   p_R.*p_L(:,l),...
            (-p_L') * FCL * p_L/2.*p_L(:,l) + FluxVal * abs(FCL)/2 * (-p_L') * (-p_L).*p_L(:,l)+...% xL
            +( p_R') * FCR * p_R/2.*p_R(:,l) + FluxVal * abs(FCR)/2 * ( p_R') *   p_R.*p_R(:,l),...
            ( p_R') * FCR * p_L/2.*p_R(:,l) + FluxVal * abs(FCR)/2 * ( p_R') * (-p_L).*p_R(:,l),...% xR
            ];
        

        %%
        % If dirichelt
        % u^-_LEFT = g(LEFT)
        % u^+_RIGHT = g(RIGHT)
        
        if strcmp(BCL,'D') %% left dirichlet
            if i==0
                TraVal = [...
                    (-p_L'*p_L) * (p_L(:,l)-p_L(:,l)),...
                    (-p_L'*p_L) * (p_L(:,l)-p_L(:,l))+...% xL
                    +( p_R') * FCR * p_R/2.*p_R(:,l) + FluxVal * abs(FCR)/2 * ( p_R') *   p_R.*p_R(:,l),...
                     ( p_R') * FCR * p_L/2.*p_R(:,l) + FluxVal * abs(FCR)/2 * ( p_R') * (-p_L).*p_R(:,l),...% xR
                    ];
            end
        end
        
        if strcmp(BCR,'D') %% right dirichlet
            if i==N-1
                TraVal = [...
                    (-p_L') * FCL * p_R/2.*p_L(:,l) + FluxVal * abs(FCL)/2 * (-p_L') *   p_R.*p_L(:,l),...
                    (-p_L') * FCL * p_L/2.*p_L(:,l) + FluxVal * abs(FCL)/2 * (-p_L') * (-p_L).*p_L(:,l)+...% xL
                    +( p_R'*p_R) * (p_R(:,l)-p_R(:,l)),...
                    ( p_R'*p_R) * (p_R(:,l)-p_R(:,l)),...% xR
                    ];
            end
        end
        

        %%
        % Adding trace value to matrix
        
        RowInd = c;
        ColInd = [c-deg,c,c+deg];
        if i == 0
            Iu = RowInd;
            Iv = ColInd(:,deg+1:end);
            Val = TraVal(:,deg+1:end);
        elseif i == N - 1
            Iu = RowInd;
            Iv = ColInd(:,1:2*deg);
            Val = TraVal(:,1:2*deg);
        else
            Iu = RowInd;
            Iv = ColInd;
            Val = TraVal;
        end
        
        tmp1D(Iu,Iv,deg*i+l) = tmp1D(Iu,Iv,deg*i+l)+Val;        
        clear Val Iu Iv
    end
end

% ======================================
% Reshape the DoF1D*DoF1D*DoF1D matrices 
%        to DoF1D^2*DoF1D
% ======================================
for i = 1:N*deg
    tmp = tmp1D(:,:,i);
    tmp = tmp';
    Mx_grad(:,i) = tmp(:);
end
% Assume matrix in y is the same as x 
My_grad = Mx_grad;

% ===== Mass Matrix ======
tmp = zeros(deg,deg,deg);
for i = 1:deg
    for j = 1:deg
        for k = 1:deg
            tmp(i,j,k) = [quad_w'* (p_val(:,i).*p_val(:,j).*p_val(:,k)) *Jacobi];
        end
    end
end

tmp1D = zeros(N*deg,N*deg,N*deg);
for i = 1:N
    ID = deg*(i-1)+[1:deg];
    tmp1D(ID,ID,ID) = tmp1D(ID,ID,ID)+tmp;
end


Mx_mass = reshape(tmp1D(:),(N*deg),(N*deg)^2);
Mx_mass = Mx_mass';
My_mass = Mx_mass;

Ind = reshape(1:(N*deg)^4,(N*deg)^2,(N*deg)^2);
count = 0;
clear ID

% Convert Mass and Grad matrices into multiwavelet basis
for i = 1:N*deg
    for j = 1:N*deg
        ix = (i-1)*N*deg;
        iy = (j-1)*N*deg;
        tmp = Ind(ix+[1:N*deg],iy+[1:N*deg]);
        ID(count+[1:(N*deg)^2],1) = tmp(:);
        count = count+(N*deg)^2;
    end
end
count = 0;

for j = 1:N*deg
    for i = 1:N*deg
        ix = (i-1)*N*deg;
        iy = (j-1)*N*deg;
        tmp = Ind(ix+[1:N*deg],iy+[1:N*deg]);
        ID_New(count+[1:(N*deg)^2],1) = tmp(:);
        count = count+(N*deg)^2;
    end
end


MM = zeros((N*deg)^2,(N*deg)^2);
for i = 1:(N*deg)
    MM([i:(N*deg):end],[i:(N*deg):end]) = FMWT_COMP;
end

ACell = repmat({FMWT_COMP}, 1, N*deg);
BigA = blkdiag(ACell{:});


% This is to mimic F*Grad*F'
MyMW_grad = MM*BigA*My_grad;
MxMW_grad = MM*BigA*Mx_grad;
MyMW_mass = MM*BigA*My_mass;
MxMW_mass = MM*BigA*Mx_mass;

% Below is to deal with some indexing issues and rearrange the index
clear ID Ind
count = 0;
Ind = reshape(1:(deg)^4,(deg)^2,(deg)^2);
for i = 1:deg
    for j = 1:deg
        ix = (i-1)*deg;
        iy = (j-1)*deg;
        tmp = Ind(iy+[1:deg],ix+[1:deg]);
        ID(count+[1:(deg)^2],1) = tmp(:);
        count = count+(deg)^2;
    end
end


% =======================
% Generate the SG Table
% =======================
Dim = 2;
[forwardHash,inverseHash] = hash_table_2D(lev,Dim,'SG');

% Assemble the SG matrix according to Hash Table
DoF_SG = numel(inverseHash);
G_sg = sparse(DoF_SG*deg^2,DoF_SG*deg^2);
for i = 1:DoF_SG
    tmp_i = inverseHash{i};

    Ix = tmp_i(end-1);
    Iy = tmp_i(end);
    
    for j = 1:DoF_SG
        tmp_j = inverseHash{j};
  
        Jx = tmp_j(end-1);
        Jy = tmp_j(end);
        
        I_ind = (2^lev)*deg*[0:deg-1]'+(Ix-1)*2^lev*deg^2+(Jx-1)*deg+[1:deg];
        J_ind = (2^lev)*deg*[0:deg-1]'+(Iy-1)*2^lev*deg^2+(Jy-1)*deg+[1:deg];
        
        val = kronmult2(MxMW_mass(I_ind(:),:),MyMW_grad(J_ind(:),:),k_vect) + ...
              kronmult2(MxMW_grad(I_ind(:),:),MyMW_mass(J_ind(:),:),h_vect);
      

        Ind_x = (i-1)*deg^2+[1:deg^2];
        Ind_y = (j-1)*deg^2+[1:deg^2];
        
        G_sg(Ind_x,Ind_y) = G_sg(Ind_x,Ind_y) + reshape(val(ID),deg^2,deg^2);
        
    end
end

% I did the truncation here!
ix = abs(G_sg)<1e-10;
G_sg(ix) = 0;
spy(G_sg)

% ==================================================
% Assemble the RHS
% Below can be optimized
% I am current doing the generation of full-grid RHS
% Then select the ones corresponding to SG
% ==================================================
[quad_x,quad_w] = lgwt(default_quad_number(deg),-1,1);
p_val = transpose( lin_legendre(quad_x,deg) * 1/sqrt(h) ); 

% for full-grid RHS calculation
count = 0;
RHS_FG = zeros((N*deg)^2,1);
for Nx =  0:N-1
    for Kx = 1:deg
        for Ny = 0:N-1
            for Ky = 1:deg
                px =  h*(quad_x/2+1/2+Nx) + xMin;
                py =  h*(quad_x/2+1/2+Ny) + xMin;
                [px,py] = meshgrid(px,py);
                val = kron(p_val(Kx,:)',p_val(Ky,:)').*kron(quad_w,quad_w).*FunF(px(:),py(:));
                RHS_FG(count+1) = sum(val)*h^2/4;
                count = count+1;
            end
        end
    end
end

RHS_FG = kronmult2(FMWT_COMP,FMWT_COMP,RHS_FG);

% Select the SG's RHS
SG_ID = zeros(size(inverseHash,2)*deg^2,1);
for i = 1:size(inverseHash,2)
    tmp_i = inverseHash{i};
    Ix = tmp_i(end-1);
    LevX = tmp_i(1);
    CelX = tmp_i(3);
 
    Iy = tmp_i(end);
    LevY = tmp_i(2);
    CelY = tmp_i(4);
    for Kx = 1:deg
        for Ky = 1:deg
            
            GlobalID = (Ix-1)*deg*(deg*N)+(Kx-1)*(deg*N)+(Iy-1)*deg+Ky;
            
            SG_ID(deg^2*(i-1)+deg*(Kx-1)+Ky) = GlobalID;
        end
    end
    
    
end
RHS_SG = RHS_FG(SG_ID);

u_sg2 = G_sg\RHS_SG;

% ==================================
% Plot u_sg on the grid
% I just copied it from ASGARD
% Only check the left bottom figures 
% ==================================
opts.deg = deg;
opts.output_grid = 'quadrature_with_end_points';
opts.max_lev = lev;
opts.lev = lev;
opts.use_oldhash = 1;
pde = TestNonSeperableCoef(opts);
for d = 1:2
    [Meval{d},nodes{d}] = matrix_plot_D(pde,opts,pde.dimensions{d});
    nodes_nodups{d} = nodes{d};
    nodes_count{d} = nodes{d}.*0+1;
end
coord = get_realspace_coords(pde,nodes);
coord_nodups = get_realspace_coords(pde,nodes_nodups);

% Plot Numerical Solution
figure
fval_realspace = wavelet_to_realspace(pde,opts,Meval,u_sg2,inverseHash);
f_realspace_nD = singleD_to_multiD(2,fval_realspace,nodes);
plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_nD);

figure;
fval_realspace = wavelet_to_realspace(pde,opts,Meval,RHS_SG,inverseHash);
f_realspace_nD = singleD_to_multiD(2,fval_realspace,nodes);
plot_fval(pde,nodes_nodups,f_realspace_nD,f_realspace_nD);


