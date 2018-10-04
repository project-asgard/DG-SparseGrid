function [gMass,Grad,Delta]=matrix_coeff_TI2(Lev,k,Lmin,Lmax,FMWT_COMP)
%=============================================================
% Generate time-independent coefficient matrices
% Vlasolv Solver:
%   Operators:  gMass: int_v v*l_i(v)*l_j(v)dv
%               GradV: int_v (l_i(v))'*l_j(v)dv
%               Grad: int_x (m_i(x))'*m_j(x)dx
% Poisson Solver:
%   Operators: Delta: int_x (m_i(x))''*m_j(x)dx
%   This equation is solved by LDG methods
% Maxwell Solver: (Missing)
%   Operators: CurlX: int_x curl(m_i(x))*m_j(x)dx
% Input: Lev, k, dim, Lmax, Vmax
% P.S. This is the full-grid version
%=============================================================
%--Quadrature
quad_num=10;
%---------------

% ------------------------------------
% compute the trace values
% note vector p_1[] is 1 by k, 
%      vector p_2[] is 1 by k
% ------------------------------------
p_1 = legendre(-1,k);
p_2 = legendre(1,k);

% ----------------------------------
% note quad_x(:) is quad_num by 1
%      quad_w(:) is quad_num by 1
%      p_val(:,:) is quad_num by k
%      Dp_val(:,:) is quad_num by k
% ----------------------------------
[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,k);
Dp_val = dlegendre(quad_x,k);

%---------------------------
% Jacobi of variable x and v
% Define Matrices
%---------------------------
nx=2^(Lev);
hx=(Lmax-Lmin)/nx;
Jacobi_x=hx;
dof_1D_x=k*nx;

nmax = 2*1024;
use_dense = (dof_1D_x <= nmax);

if (use_dense),
  Grad=zeros(dof_1D_x,dof_1D_x);
  size_Delta = 2*dof_1D_x;
  Delta=zeros(size_Delta,size_Delta);
else
  Grad=sparse(dof_1D_x,dof_1D_x);
  Delta=sparse(2*dof_1D_x,2*dof_1D_x);
end;


% nx=2^(Lev_v);
% hx=2*Vmax/nx;
% Jacobi_x=hx;
% dof_1D_x=k*nx;
% 
if (use_dense),
  gMass=zeros(dof_1D_x,dof_1D_x);
%   GradV=zeros(dof_1D_x,dof_1D_x);
else
  gMass=sparse(dof_1D_x,dof_1D_x);
%   GradV=sparse(dof_1D_x,dof_1D_x);
end;

%======================================
% Matrices related to x variable
% Grad and Delta
%======================================
% compute (u',v)+1/2*u^{-}[v^{+}-v^{-}]
% generate 1D matrix for DG
for Lx=0:nx-1

    %---------------------------------------------
    % Matrix Grad and EMassX
    %---------------------------------------------

    % -----------------------
    % note val(:,:) is k by k
    % -----------------------
    val=(1/hx)*[Dp_val'*(quad_w.*p_val)];

    i1 = k*Lx+1;
    i2 = k*(Lx+1);
    % -----------------------------
    % note  (i2-i1+1) is equal to k
    % -----------------------------
    

    if (use_dense),
      Grad(i1:i2,i1:i2) = Grad(i1:i2,i1:i2) + val(1:k,1:k); 
    else
    % ----------------------------
    % note Iv(:,:) is  k by k
    % Iv = [  i1, i1+1, ..., i2; 
    %         i1, i1+1, ..., i2;
    %         ...
    %         i1, i1+1, ..., i2]
    %
    % Iu = [ i1,   i1,   ..., i1;
    %        i1+1, i1+1, ..., i1+1;
    %        ...
    %        i2,   i2,   ..., i2]
    % ----------------------------
      Iv = meshgrid(i1:i2);
      Iu = transpose(Iv);
      Grad=Grad+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
    end;
    
    if (use_dense),
      % --------------------------------------------------------------------------
      % Delta = Delta + sparse( Iu, dof_1D_x+Iv, val, size_Delta, size_Delta);
      % --------------------------------------------------------------------------
      Delta(i1:i2, (i1+dof_1D_x):(i2+dof_1D_x)) = ...
            Delta(i1:i2,(i1+dof_1D_x):(i2+dof_1D_x)) + val(1:k,1:k);

      % --------------------------------------------------------------------------
      % Delta = Delta + sparse( dof_1D_x+Iu, Iv, val, size_Delta, size_Delta);
      % --------------------------------------------------------------------------
      Delta( (i1+dof_1D_x):(i2+dof_1D_x), i1:i2 ) = ...
            Delta( (i1+dof_1D_x):(i2+dof_1D_x), i1:i2) + val(1:k,1:k);

      % ----------------------------------------------------------------------
      % Delta = Delta + sparse( Iu, Iv, eye(k,k), size_Delta, size_Delta);
      % ----------------------------------------------------------------------
      Delta(i1:i2,i1:i2) = Delta(i1:i2,i1:i2) + eye(k,k);
      
    else
      Delta=Delta+sparse([Iu,dof_1D_x+Iu,Iu],[dof_1D_x+Iv,Iv,Iv],...
        [val,val,diag(ones(1,k))],2*dof_1D_x,2*dof_1D_x);
    end;
    
    
    % ------------------
    % c=k*Lx+1:k*(Lx+1);
    % ------------------
    c1 = k*Lx+1;
    c2 = k*(Lx+1);
    c = c1:c2;

    % ------------------
    % p=k*(Lx-1)+1:k*Lx;
    % ------------------
    p1 = k*(Lx-1)+1;
    p2 = k*Lx;
    p = p1:p2;

    % ----------------------
    % l=k*(Lx+1)+1:k*(Lx+2);
    % ----------------------
    l1 = k*(Lx+1)+1;
    l2 = k*(Lx+2);
    l = l1:l2;
    
    val=1/hx*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
        p_2'*p_2/2   p_2'*p_1/2];     % for x2
    
    val_u=1/hx*[-p_1'*p_1, p_2'*p_1];
    val_s=1/hx*[-p_1'*p_2, p_2'*p_2];
    
    
    if (use_dense),
    % -----------------------------------------------
    % setup ranges for row indices and column indices
    % -----------------------------------------------
     istart = zeros(4,1);
     iend = zeros(4,1);
     jstart = zeros(4,1);
     jend = zeros(4,1);


     % ----------------------------------------------------------
     % Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
     % ----------------------------------------------------------
   
     istart(1) = c1; iend(1) = c2;
     istart(2) = c1; iend(2) = c2;
     istart(3) = c1; iend(3) = c2;
     istart(4) = c1; iend(4) = c2;

    else
      Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
    end;


    
    if Lx<nx-1 && Lx>0
        
     if (use_dense),
        % ----------------------------------------------------
        % Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
        % ----------------------------------------------------
       
        jstart(1) = p1; jend(1) = p2;
        jstart(2) = c1; jend(2) = c2;
        jstart(3) = c1; jend(3) = c2;
        jstart(4) = l1; jend(4) = l2;
      else
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
      end;
    
    elseif Lx==0
        
      if (use_dense),
        % ----------------------------------------------------------------------
        % Iu=[meshgrid([k*(nx-1)+1:k*(nx)]),meshgrid(c),meshgrid(c),meshgrid(l)];
        % ----------------------------------------------------------------------
        jstart(1) = (k*(nx-1)+1); jend(1) = k*(nx);
        jstart(2) = c1;           jend(2) = c2;
        jstart(3) = c1;           jend(3) = c2;
        jstart(4) = l1;           jend(4) = l2;

      else
        Iu=[meshgrid([k*(nx-1)+1:k*(nx)]),meshgrid(c),meshgrid(c),meshgrid(l)];
      end;
   
    elseif Lx==nx-1
        
      if (use_dense),
        % ---------------------------------------------------------
        % Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
        % ---------------------------------------------------------
        jstart(1) = p1; jend(1) = p2;
        jstart(2) = c1; jend(2) = c2;
        jstart(3) = c1; jend(3) = c2;
        jstart(4) =  1; jend(4) = k;

      else
        Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
      end;
    
    end
    
    if (use_dense),
      % ---------------------
      % note subtraction used
      % ---------------------
      ia1 = istart(1); 
      ia2 = iend(1);
      ja1 = jstart(1); 
      ja2 = jend(1);

      Grad(ia1:ia2, ja1:ja2) = Grad(ia1:ia2, ja1:ja2) - ...
             val( 1:k, 1:k);

      ia1 = istart(2); 
      ia2 = iend(2);
      ja1 = jstart(2); 
      ja2 = jend(2);

      Grad(ia1:ia2, ja1:ja2) = Grad(ia1:ia2, ja1:ja2) - ...
             val( 1:k, k+(1:k));


      ia1 = istart(3); 
      ia2 = iend(3);
      ja1 = jstart(3); 
      ja2 = jend(3);
      
      Grad(ia1:ia2, ja1:ja2) = Grad(ia1:ia2, ja1:ja2) - ...
             val( 1:k, 2*k+(1:k));

      ia1 = istart(4); 
      ia2 = iend(4);
      ja1 = jstart(4); 
      ja2 = jend(4);

      Grad(ia1:ia2, ja1:ja2) = Grad(ia1:ia2, ja1:ja2) - ...
             val( 1:k, 3*k+(1:k));

    else
      Grad=Grad-sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
    end;

    if (use_dense),
% ------------------------------------------
%      Delta=Delta - ...
%       sparse([Iv(:,1:2*k),Iv(:,1:2*k)+dof_1D_x],...
%              [Iu(:,2*k+1:end)+dof_1D_x,Iu(:,1:2*k)],...
%              [val_u,val_s],2*dof_1D_x,2*dof_1D_x);
% ------------------------------------------


%        ---------------------------------------
%        Delta = Delta - ...
%          sparse( Iv(1:k, 1:k), ...
%                  Iu(1:k, 2*k + (1:k)) + dof_1D_x, ...
%                  val_u(1:k,1:k), ...
%                  size_Delta,size_Delta);
%        ---------------------------------------
         ia1 = istart(1);
         ia2 = iend(1);
         ja1 = jstart(3) + dof_1D_x;
         ja2 = jend(3)   + dof_1D_x;

         Delta(ia1:ia2, ja1:ja2) = Delta(ia1:ia2, ja1:ja2) - ...
              val_u(1:k,1:k);

%        ---------------------------------------
%        Delta = Delta - ...
%          sparse( Iv(1:k, k + (1:k)), ...
%                  Iu(1:k, 3*k  (1:k)) + dof_1D_x, ...
%                  val_u(1:k, k+(1:k)), ...
%                  size_Delta,size_Delta);
%        ---------------------------------------
         ia1 = istart(2);
         ia2 = iend(2);
         ja1 = jstart(4) + dof_1D_x;
         ja2 = jend(4)   + dof_1D_x;
  
         Delta(ia1:ia2,ja1:ja2) = Delta(ia1:ia2,ja1:ja2) - ...
              val_u(1:k, k + (1:k) );

%        ---------------------------------------
%        Delta = Delta - ...
%          sparse( Iv(1:k, 1:k) + dof_1D_x, ...
%                  Iu(1:k, 1:k), ...
%                  val_s(1:k, 1:k), ...
%                  size_Delta,size_Delta);
%        ---------------------------------------
         ia1 = istart(1) + dof_1D_x;
         ia2 = iend(1)   + dof_1D_x;
         ja1 = jstart(1);
         ja2 = jend(1);

         Delta(ia1:ia2, ja1:ja2) = Delta(ia1:ia2, ja1:ja2) - ...
              val_s(1:k, 1:k);

%        ---------------------------------------
%        Delta = Delta - ...
%          sparse( Iv(1:k, k + (1:k)) + dof_1D_x, ...
%                  Iu(1:k, k + (1:k)), ...
%                  val_s(1:k, k+(1:k)), ...
%                  size_Delta,size_Delta);
%        ---------------------------------------
         ia1 = istart(2) + dof_1D_x; 
         ia2 = iend(2)   + dof_1D_x;
         ja1 = jstart(2);
         ja2 = jend(2);
        
         Delta(ia1:ia2, ja1:ja2) = Delta(ia1:ia2,ja1:ja2) - ...
              val_s(1:k, k+(1:k) );

    else
      Delta=Delta+sparse([Iv(:,1:2*k),Iv(:,1:2*k)+dof_1D_x],...
        ...[Iu(:,1:2*k)+dof_1D_x,Iu(:,2*k+1:end)],...
        [Iu(:,2*k+1:end)+dof_1D_x,Iu(:,1:2*k)],...
        -[val_u,val_s],2*dof_1D_x,2*dof_1D_x);
    end;
    
    
end

% Handle B.C. for Poisson solver
Delta(dof_1D_x+1,:)=0;
Delta(dof_1D_x+1,dof_1D_x+[1:k])=sqrt(1/hx)*legendre(-1,k);

% % iend = size_Delta;
% % Delta(iend,1:size_Delta)=0;
% % % ---------------------------------------------------
% % % Delta(iend,iend-k+[1:k])=sqrt(1/hx)*legendre(1,k);
% % % ---------------------------------------------------
% % leg_1_k = zeros(1,k);
% % leg_1_k(1:k) = legendre(1,k);
% % k1 = (iend-k+1);
% % k2 = iend-k+k;
% % Delta(iend,k1:k2) = sqrt(1/hx)*leg_1_k(1:k);

%======================================
% Matrices related to v variable
% gMass and GradV
%======================================
% (vf(v),w(v))_Kv
for Lv=0:nx-1
    xi_v=(( (quad_x+1)/2+Lv)*hx-Lmin); % mapping from [-1,1] to physical domain
    %---------------------------------------------
    % Matrix
    %---------------------------------------------
    % value of local matrix
    val_loc=p_val'*(p_val.*xi_v.*quad_w)*Jacobi_x/2/hx;

    if (use_dense),
      i1 = k*Lv+1;
      i2 = k*(Lv+1);
      % --------------------
      % note (i2-i1+1) == k
      % --------------------
      gMass(i1:i2, i1:i2) = gMass(i1:i2, i1:i2) + ...
            val_loc(1:k,1:k);
    else
      Iu=meshgrid(k*Lv+1:k*(Lv+1));
      gMass=gMass+sparse(Iu',Iu,val_loc,dof_1D_x,dof_1D_x);
    end;
    
%     val=1/hx*[Dp_val'*(quad_w.*p_val)];
%     if (use_dense),
%       i1 = k*Lv+1;
%       i2 = k*(Lv+1);
%       % -------------------
%       % note (i2-i1+1) == k
%       % -------------------
%       GradV(i1:i2, i1:i2) = GradV(i1:i2, i1:i2) +  ...
%            val(1:k,1:k);
%     else
%       Iu=[meshgrid(k*Lv+1:k*(Lv+1))]';
%       Iv=[meshgrid(k*Lv+1:k*(Lv+1))];
%       GradV=GradV+sparse(Iu,Iv,val,dof_1D_x,dof_1D_x);
%     end;
    
%     % -----------------
%     % c=k*Lv+1:k*(Lv+1);
%     % -----------------
%     c1 = k*Lv+1;
%     c2 = k*(Lv+1);
%     c = c1:c2;
% 
% 
%     % ------------------
%     % p=k*(Lv-1)+1:k*Lv;
%     % ------------------
%     p1 = k*(Lv-1)+1;
%     p2 = k*Lv;
%     p = p1:p2;
% 
%     % ----------------------
%     % l=k*(Lv+1)+1:k*(Lv+2);
%     % ----------------------
%     l1 = k*(Lv+1)+1;
%     l2 = k*(Lv+2);
%     l = l1:l2;
%     
%     val=(1/hx)*[-p_1'*p_2/2  -p_1'*p_1/2,...   % for x1
%         p_2'*p_2/2   p_2'*p_1/2];     % for x2
%     
%     istart = zeros(4,1);
%     iend = zeros(4,1);
%     jstart = zeros(4,1);
%     jend = zeros(4,1);
% 
%     if Lv<nx-1 && Lv>0
%      if (use_dense)
%         % -----------------------------------------------------
%         % Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
%         % -----------------------------------------------------
%         jstart(1) = p1; jend(1) = p2;
%         jstart(2) = c1; jend(2) = c2;
%         jstart(3) = c1; jend(3) = c2;
%         jstart(4) = l1; jend(4) = l2;
% 
%         % --------------------------------------------------------
%         % Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
%         % --------------------------------------------------------
%         istart(1) = c1; iend(1) = c2;
%         istart(2) = c1; iend(2) = c2;
%         istart(3) = c1; iend(3) = c2;
%         istart(4) = c1; iend(4) = c2;
% 
% 
%      else   
%         Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid(l)];
%         Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
%      end;
%         
%     elseif Lv==0
%       if (use_dense),
%         % ----------------------------------------------------------------------
%         % Iu=[meshgrid([k*(nx-1)+1:k*(nx)]),meshgrid(c),meshgrid(c),meshgrid(l)];
%         % ----------------------------------------------------------------------
%         jstart(1) = k*(nx-1)+1; jend(1) = k*(nx);
%         jstart(2) = c1;         jend(2) = c2;
%         jstart(3) = c1;         jend(3) = c2;
%         jstart(4) = l1;         jend(4) = l2;
%         
%         % ----------------------------------------------------------------------
%         % Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
%         % ----------------------------------------------------------------------
%         istart(1) = c1;  iend(1) = c2;
%         istart(2) = c1;  iend(2) = c2;
%         istart(3) = c1;  iend(3) = c2;
%         istart(4) = c1;  iend(4) = c2;
% 
%       else 
%         Iu=[meshgrid([k*(nx-1)+1:k*(nx)]),meshgrid(c),meshgrid(c),meshgrid(l)];
%         Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
%       end;
%         
%     elseif Lv==nx-1
%         
%      if (use_dense),
%         % ---------------------------------------------------------
%         % Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
%         % ---------------------------------------------------------
%         jstart(1) = p1; jend(1) = p2;
%         jstart(2) = c1; jend(2) = c2;
%         jstart(3) = c1; jend(3) = c2;
%         jstart(4) = 1;  jend(4) = k;
% 
%         % ---------------------------------------------------------
%         % Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
%         % ---------------------------------------------------------
%         istart(1) = c1; iend(1) = c2;
%         istart(2) = c1; iend(2) = c2;
%         istart(3) = c1; iend(3) = c2;
%         istart(4) = c1; iend(4) = c2;
%      else
%         Iu=[meshgrid(p),meshgrid(c),meshgrid(c),meshgrid([1:k])];
%         Iv=[meshgrid(c)',meshgrid(c)',meshgrid(c)',meshgrid(c)'];
%      end;
%     end
% 
%     if (use_dense),
%        % ---------------------
%        % note subtraction used
%        % ---------------------
%        ia1 = istart(1); ia2 = iend(1);
%        ja1 = jstart(1); ja2 = jend(1);
% 
%        GradV(ia1:ia2, ja1:ja2) = GradV(ia1:ia2,ja1:ja2) - ...
%             val(1:k, 1:k);
% 
%        ia1 = istart(2); ia2 = iend(2);
%        ja1 = jstart(2); ja2 = jend(2);
%        GradV(ia1:ia2, ja1:ja2) = GradV(ia1:ia2,ja1:ja2) - ...
%             val(1:k, k + (1:k));
% 
%        ia1 = istart(3); ia2 = iend(3);
%        ja1 = jstart(3); ja2 = jend(3);
%        GradV(ia1:ia2, ja1:ja2) = GradV(ia1:ia2,ja1:ja2) - ...
%             val(1:k, 2*k + (1:k));
% 
%        ia1 = istart(4); ia2 = iend(4);
%        ja1 = jstart(4); ja2 = jend(4);
%        GradV(ia1:ia2, ja1:ja2) = GradV(ia1:ia2,ja1:ja2) - ...
%             val(1:k, 3*k + (1:k));
%     else
%        GradV=GradV-sparse(Iv,Iu,val,dof_1D_x,dof_1D_x);
%     end;
end

%***************************************
% Following is for Multiwavelet DG Basis
% Simple way to conxert??
%***************************************
% Transfer to multi-DG bases
gMass = FMWT_COMP*gMass*FMWT_COMP';
Grad = FMWT_COMP*Grad*FMWT_COMP';
% GradV = FMWT_COMP_v*GradV*FMWT_COMP_v';

use_blkdiag = 0;
if (use_blkdiag),
  Delta = blkdiag(FMWT_COMP,FMWT_COMP)*...
                Delta*...
           blkdiag(FMWT_COMP',FMWT_COMP');
else
 % ----------------------------------------
 % note blkdiag(A,B) creates a block diagonal matrix
 % [A, 0;
 % [0, B]
 %
 % let F = FMWT_COMP_x
 %
 % Delta = [F, 0;   [D11,  D12;   [F', 0;
 %           0, F] * [D21,  D22] *  0 , F']
 % Delta = [ F*D11*F',    F*D12*F';
 %            F*D21*F',    F*D22*F']
 % computed as
 % step 1:  Delta = blkdiag(F,F) * Delta
 % to form  [F * D11, F * D12;
 %           F * D21, F * D22 ]
 %
 % step 2:  Delta = Delta * blkdiag(F',F')
 % to form [ (F*D11)*F',   (F*D12)*F';
 %           (F*D21)*F',   (F*D22)*F']
 % ----------------------------------------

 % ---------------------------------
 % note shape of Delta matrix 
 % is  (2*dof_1D_x)  by (2*dof_1D_x)
 % ---------------------------------
 m = dof_1D_x;
 F = FMWT_COMP;

 % --------------------------------
 % step 1: multiply by blkdiag(F,F) on left
 % --------------------------------
 
 Delta(1:m, 1:(2*m)) = F * Delta(1:m,1:(2*m));
 Delta((m+1):(2*m), 1:(2*m)) = F * Delta( (m+1):(2*m), 1:(2*m) );

 % -----------------------------------
 % multiply by blkdiag(F',F') on right
 % -----------------------------------

 Delta(1:(2*m), 1:m) = Delta(1:(2*m),1:m) * F';
 Delta(1:(2*m), (m+1):(2*m)) = Delta(1:(2*m),(m+1):(2*m)) * F';
end;
