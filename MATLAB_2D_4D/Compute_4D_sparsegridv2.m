% Computing on the sparse grid for Poisson Eq
% 4-dimensional calculation
%------------------------------------------------
% A_s: Matrix
% b_s: RHS
% sol_s: Solution
% uu_s: Interpolation
% Loop sum_level
%------------------------------------------------

dim=4;

% Key1DMesh
% ['Generate 1D keys']
% tic
nx=[];px=[];kx=[];
for Lx=0:n
    for Px=0:2^max(0,Lx-1)-1
        for Kx=1:k
            nx=[nx;Lx];
            px=[px;Px];
            kx=[kx;Kx];
            
        end
    end
end
%
Key1dMesh=[nx,px,kx];
% toc

['Generate ',num2str(dim),'-D Hash Table']
run hashTest4D.m


A_s = sparse(dof_sparse,dof_sparse);
b_s = sparse(dof_sparse,1);
sol_s = sparse(dof_sparse,1);
uu_s=sparse(dof_sparse,1);

count = 1;
kron_flops = 0;
kron_nnz = 0;


% Method
count = 1;
clear A_encode;
for sum_level=0:n
    for i1_level=0:sum_level
        for i2_level=0:sum_level-i1_level
            for i3_level=0:sum_level-i1_level-i2_level
                i4_level=sum_level-i1_level-i2_level-i3_level;
                
                I1=Index_1D(k,i1_level);
                I2=Index_1D(k,i2_level);
                I3=Index_1D(k,i3_level);
                I4=Index_1D(k,i4_level);
                
                
                
                key_i=GenerateKey4Dv2(I1,I2,I3,I4,Key1dMesh);

                clear Index_I;
                for iii=1:size(key_i,1)
                    Index_I(iii,1)=database.(sprintf('i%g_',key_i(iii,:)));
                end
                
                
                
                % Term 1: S*I*I*I+I*S*I*I+I*I*S*I+I*I*I*S at the diagonal entry
                tmp=kron(kron(kron(Stiff_1D(I1,I1),M_mass(I2,I2)),M_mass(I3,I3)),M_mass(I4,I4))+...
                    kron(kron(kron(M_mass(I1,I1),Stiff_1D(I2,I2)),M_mass(I3,I3)),M_mass(I4,I4))+...
                    kron(kron(kron(M_mass(I1,I1),M_mass(I2,I2)),Stiff_1D(I3,I3)),M_mass(I4,I4))+...
                    kron(kron(kron(M_mass(I1,I1),M_mass(I2,I2)),M_mass(I3,I3)),Stiff_1D(I4,I4));
                

                [II,JJ]=meshgrid(Index_I,Index_I);
                A_s=A_s+sparse(II,JJ,tmp,dof_sparse,dof_sparse);
                
                tmp=kron(kron(kron(b(I1),b(I2)),b(I3)),b(I4));
                b_s(Index_I)=b_s(Index_I)+tmp;
                uu_s(Index_I)=uu_s(Index_I)+...
                    kron(kron(kron(coef_MW(I1),coef_MW(I2)),coef_MW(I3)),coef_MW(I4));
                
                % save matrices to files
                A_encode{count}.A1=Stiff_1D(I1,I1);
                A_encode{count}.A2=M_mass(I2,I2);
                A_encode{count}.A3=M_mass(I3,I3);
                A_encode{count}.A4=M_mass(I4,I4);

                A_encode{count}.IndexI = Index_I;
                A_encode{count}.IndexJ = Index_I;

                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                count = count + 1;

                
                A_encode{count}.A1=M_mass(I1,I1);
                A_encode{count}.A2=Stiff_1D(I2,I2);
                A_encode{count}.A3=M_mass(I3,I3);
                A_encode{count}.A4=M_mass(I4,I4);

                A_encode{count}.IndexI = Index_I;
                A_encode{count}.IndexJ = Index_I;

                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                count = count + 1;

                
                A_encode{count}.A1=M_mass(I1,I1);
                A_encode{count}.A2=M_mass(I2,I2);
                A_encode{count}.A3=Stiff_1D(I3,I3);
                A_encode{count}.A4=M_mass(I4,I4);

                A_encode{count}.IndexI = Index_I;
                A_encode{count}.IndexJ = Index_I;


                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );

                kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                count = count + 1;


                
                A_encode{count}.A1=M_mass(I1,I1);
                A_encode{count}.A2=M_mass(I2,I2);
                A_encode{count}.A3=M_mass(I3,I3);
                A_encode{count}.A4=Stiff_1D(I4,I4);

                A_encode{count}.IndexI=Index_I;
                A_encode{count}.IndexJ=Index_I;

                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );


                kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                
                count = count + 1;


                
                % Term 2: S*I*I*I--Assume j2_level==i2_level
                % j3_level==i3_level j4_level==i4_level
                j2_level=i2_level;j3_level=i3_level;j4_level=i4_level;
                for j1_level=0:i1_level-1
                    
                    J1=Index_1D(k,j1_level);
                    J2=I2;J3=I3;J4=I4;
                    

                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);

                    clear Index_J;
                    for jjj=1:size(key_j,1)
                       Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
                    end
                    
                    
                    [II,JJ]=meshgrid(Index_J,Index_I);
                    
                    tmp=kron(kron(kron(Stiff_1D(J1,I1),M_mass(J2,I2)),M_mass(J3,I3)),M_mass(J4,I4));
                    
                    A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=Stiff_1D(J1,I1);
                    A_encode{count}.A2=M_mass(J2,I2);
                    A_encode{count}.A3=M_mass(J3,I3);
                    A_encode{count}.A4=M_mass(J4,I4);

                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;
                    
                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    
                    


                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    count = count + 1;


                    
                    A_encode{count}.A1=Stiff_1D(J1,I1)';
                    A_encode{count}.A2=M_mass(J2,I2)';
                    A_encode{count}.A3=M_mass(J3,I3)';
                    A_encode{count}.A4=M_mass(J4,I4)';

                    A_encode{count}.IndexI=Index_I;
                    A_encode{count}.IndexJ=Index_J;

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    
                    

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    
                    count=count+1;
                end
                
                % Term 3: I*S*I*I--Assume j1_level==i1_level
                % j3_level==i3_level j4_level==i4_level
                j1_level=i1_level;j3_level=i3_level;j4_level=i4_level;
                for j2_level=0:i2_level-1
                    
                    J2=Index_1D(k,j2_level);
                    J1=I1;J3=I3;J4=I4;
                    
                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);

                    clear Index_J;
                    for jjj=1:size(key_j,1)
                        Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
                    end
                    
                    
                    
                    [II,JJ]=meshgrid(Index_J,Index_I);
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),Stiff_1D(J2,I2)),M_mass(J3,I3)),M_mass(J4,I4));
                    
                    A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    
                    A_encode{count}.A1=M_mass(J1,I1);
                    A_encode{count}.A2=Stiff_1D(J2,I2);
                    A_encode{count}.A3=M_mass(J3,I3);
                    A_encode{count}.A4=M_mass(J4,I4);

                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    count = count + 1;
                    
                    
                    A_encode{count}.A1=M_mass(J1,I1)';
                    A_encode{count}.A2=Stiff_1D(J2,I2)';
                    A_encode{count}.A3=M_mass(J3,I3)';
                    A_encode{count}.A4=M_mass(J4,I4)';

                    A_encode{count}.IndexI=Index_I;
                    A_encode{count}.IndexJ=Index_J;

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );


                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);

                    count=count+1;
                end
                
                % Term 4: I*I*S*I--Assume j1_level==i1_level
                % j2_level==i2_level j4_level==i4_level
                j1_level=i1_level;j2_level=i2_level;j4_level=i4_level;
                for j3_level=0:i3_level-1
                    
                    J3=Index_1D(k,j3_level);
                    J1=I1;J2=I2;J4=I4;
                    
                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);

                    clear Index_J;
                    for jjj=1:size(key_j,1)
                        Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
                    end
                    
                    
                    
                    [II,JJ]=meshgrid(Index_J,Index_I);
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),M_mass(J2,I2)),Stiff_1D(J3,I3)),M_mass(J4,I4));
                    A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    
                    A_encode{count}.A1=M_mass(J1,I1);
                    A_encode{count}.A2=M_mass(J2,I2);
                    A_encode{count}.A3=Stiff_1D(J3,I3);
                    A_encode{count}.A4=M_mass(J4,I4);

                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    count = count + 1;
                    
                    
                    A_encode{count}.A1=M_mass(J1,I1)';
                    A_encode{count}.A2=M_mass(J2,I2)';
                    A_encode{count}.A3=Stiff_1D(J3,I3)';
                    A_encode{count}.A4=M_mass(J4,I4)';

                    A_encode{count}.IndexI=Index_I;
                    A_encode{count}.IndexJ=Index_J;

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    
                    count=count+1;
                end
                
                % Term 5: I*I*I*S--Assume j1_level==i1_level
                % j2_level==i2_level j3_level==i3_level
                j1_level=i1_level;j2_level=i2_level;j3_level=i3_level;
                for j4_level=0:i4_level-1
                    
                    J4=Index_1D(k,j4_level);
                    J1=I1;J2=I2;J3=I3;
                    

                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);

                    clear Index_J;
                    for jjj=1:size(key_j,1)
                        Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
                    end
 
                    [II,JJ]=meshgrid(Index_J,Index_I);
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),M_mass(J2,I2)),M_mass(J3,I3)),Stiff_1D(J4,I4));
                    A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    
                    A_encode{count}.A1=M_mass(J1,I1);
                    A_encode{count}.A2=M_mass(J2,I2);
                    A_encode{count}.A3=M_mass(J3,I3);
                    A_encode{count}.A4=Stiff_1D(J4,I4);

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    
                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    count = count + 1;

                    
                    
                    A_encode{count}.A1=M_mass(J1,I1)';
                    A_encode{count}.A2=M_mass(J2,I2)';
                    A_encode{count}.A3=M_mass(J3,I3)';
                    A_encode{count}.A4=Stiff_1D(J4,I4)';

                    A_encode{count}.IndexI=Index_I;
                    A_encode{count}.IndexJ=Index_J;

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    
                    count=count+1;
                end
                
            end
        end
    end
    
end
clear Index_I;
clear Index_J;

save(['./Data/A_4D_encode.mat'],'A_encode');

% -----------------------
% ensure A_s is symmetric
% -----------------------
A_s  = (A_s + A_s')/2;

disp(sprintf('Np=%d,k=%d, n=%d, kron_flops=%g, kron_nnz=%g, nnz(A_s)=%g ', ...
              Np, k, size(A_s,1), ...
              kron_flops,    kron_nnz, nnz(A_s)  ));

use_dense_matrix = 1;
if (use_dense_matrix),
  ratio = 0.25;
  for count=1:length(A_encode),

     is_very_sparse = (nnz(A_encode{count}.A1) < ...
                      ratio*prod(size(A_encode{count}.A1)));
     if (~is_very_sparse),
       A_encode{count}.A1 = full(A_encode{count}.A1);
     end;

     is_very_sparse = (nnz(A_encode{count}.A2) < ...
                      ratio*prod(size(A_encode{count}.A2)));
     if (~is_very_sparse),
       A_encode{count}.A2 = full(A_encode{count}.A2);
     end;

     is_very_sparse = (nnz(A_encode{count}.A3) < ...
                      ratio*prod(size(A_encode{count}.A3)));
     if (~is_very_sparse),
       A_encode{count}.A3 = full(A_encode{count}.A3);
     end;

     is_very_sparse = (nnz(A_encode{count}.A4) < ...
                      ratio*prod(size(A_encode{count}.A4)));
     if (~is_very_sparse),
       A_encode{count}.A4 = full(A_encode{count}.A4);
     end;

  end;
end;

% -----------------------------------
% estimate work using fixed kron strategy and 
% strategy to minimize work
% assume all matrices A1..A4 are dense
% -----------------------------------
flops_fixed = 0.0;
flops_min = 0.0;
for icount=1:length(A_encode),
    nrow1 = size(A_encode{icount}.A1,1);
    nrow2 = size(A_encode{icount}.A2,1);
    nrow3 = size(A_encode{icount}.A3,1);
    nrow4 = size(A_encode{icount}.A4,1);
    
    ncol1 = size(A_encode{icount}.A1,2);
    ncol2 = size(A_encode{icount}.A2,2);
    ncol3 = size(A_encode{icount}.A3,2);
    ncol4 = size(A_encode{icount}.A4,2);

    rc = [nrow1, ncol1; ...
          nrow2, ncol2; ...
          nrow3, ncol3; ...
          nrow4, ncol4 ];
    rc = transpose(rc);

    flops1 = kron_cost_fixed(rc);
    [flops2,isplit,imethod] = kron_minflops( rc );

    flops_fixed = flops_fixed + flops1;
    flops_min = flops_min + flops2;
end;
disp(sprintf('flops for fixed kron strategy is %g', flops_fixed));
disp(sprintf('flops for min kron strategy is %g', flops_min));
    
    
  

figure;
spy(A_s)
title(sprintf('4D problem, Np=%d,k=%d, n=%d,nnz=%g',...
    Np, k, ...
    size(A_s,1),nnz(A_s)  ));

% -------------------------
% estimate condition number
% -------------------------
%eigs(A_s,2,'BE')

function y = apply_As(x, varargin)
  As = varargin{1};
  y = As * x;
endfunction

function y = apply_Aencode(x, varargin)
  
  Aencode = varargin{1};
  N = prod(size(Aencode));
  y = zeros(size(x));
  for count = 1: N
    
    Index_I = Aencode{count}.IndexI;
    Index_J = Aencode{count}.IndexJ;

    A1 = Aencode{count}.A1;
    A2 = Aencode{count}.A2;
    A3 = Aencode{count}.A3;
    A4 = Aencode{count}.A4;


    y(Index_I) = y(Index_I) + kron_mult4(A1,A2,A3,A4, x(Index_J));
  end;

endfunction

testX = 2*rand(size(b_s))-1;
tic
Y1 = apply_As(testX, A_s);
disp(sprintf('time for apply_As %g ', toc ));

tic
Y2 = apply_Aencode(testX, A_encode);
disp(sprintf('time for apply_Aencode %g', toc));

disp(sprintf('norm(Y1-Y2,1)=%g', norm(Y1-Y2,1)));
disp(sprintf('norm(Y1,1)=%g', norm(Y1,1) ));
disp(sprintf('relative error is %g', ...
    norm(Y1-Y2,1)/norm(Y1,1)  ));


% tic
is_small_A_s = (size(A_s,1) <= 4*1024);
use_direct_solve = is_small_A_s;
if (use_direct_solve),
  sol_s = A_s\b_s*pi^2*4;
else
  x0 = 0*b_s;
  tol = 1e-9;
  maxit = min( 1000, size(A_s,1)+10);

  use_ilu = 0;
  if (use_ilu),
    % ---------------------------------
    % use Incomplete LU  factorization as
    % preconditioner
    % ---------------------------------
    time_ilu = -time();
    ilu_opts.type = 'nofill';
    % ----------------
    % M1 is L, M2 is U
    % ----------------
    [M1,M2] = ilu(A_s);   
    time_ilu = time_ilu + time();
    disp(sprintf('time for ILU=%g sec ', time_ilu ));

   else
    % ---------------------------
    % use simple diagonal scaling
    % ---------------------------
    n = size(A_s,1);
    d = diag(A_s,0);
    d( abs(d) < eps ) = 1;
    M1 = sparse( 1:n,1:n, d,  n,n);
    M2 = [];
   end;

  use_function_handle = true
  use_Aencode = false
  time_iter = -time();
  if (use_function_handle),
    if (use_Aencode),
      myfun = "apply_Aencode";
      myarg = A_encode;
    else
      myfun = "apply_As";
      myarg = A_s;
    end;
    [sol_s,flag,relres,iter,resvec,eigest] = pcg( ...
       myfun, b_s*pi^2*4, tol, maxit, M1, M2, x0, myarg );
  else
    [sol_s,flag,relres,iter,resvec,eigest] = pcg( ...
        A_s, b_s*pi^2*4, tol, maxit, M1, M2, x0 );
  end;
  time_iter = time_iter + time();
  disp(sprintf('time for iterative method=%g ', time_iter));
  isok = (flag == 0);
  if (isok),
    disp(sprintf('convergence after %d iterations',iter));
  else
    disp(sprintf('NO convergence after %d iterations',iter));
  end;

  condest = eigest(2) / eigest(1);
  disp(sprintf('eigest=(%g, %g), condest=%g', ...
          eigest(1),eigest(2), condest ));


  figure;
  semilogy( resvec  );
  title('converge of iterative method ');
end;
  
% toc

save(['./Data/A_4D_encode.mat'],'A_encode');


['Done of Solution']
% check error
norm(sol_s-uu_s)

function Ix=Index_1D(k,level)

if level==0
    Ix=[1:k];
else
    Ix=[k*2^(level-1)+1:k*2^level];
end

end

function key=GenerateKey4Dv2_org(I1,I2,I3,I4,Key1dMesh)

tmp_1=Key1dMesh(I1,:);
tmp_2=Key1dMesh(I2,:);
tmp_3=Key1dMesh(I3,:);
tmp_4=Key1dMesh(I4,:);

n1=size(tmp_1,1);
n2=size(tmp_2,1);
n3=size(tmp_3,1);
n4=size(tmp_4,1);
dim=4;
ntol=n1*n2*n3*n4;
key=zeros(ntol,3*dim);
count=1;
for i1=1:n1
    for i2=1:n2
        for i3=1:n3
            for i4=1:n4
                key(count,[1:dim:3*dim])=tmp_1(i1,:);
                key(count,[2:dim:3*dim])=tmp_2(i2,:);
                key(count,[3:dim:3*dim])=tmp_3(i3,:);
                key(count,[4:dim:3*dim])=tmp_4(i4,:);
                count=count+1;
            end
        end
    end
end

key(:,3*dim-dim+1:3*dim)=key(:,3*dim-dim+1:3*dim)-1;
end

function key=GenerateKeyd(Id,Key1dMesh)

dim=size(Id,1);
n=zeros(dim,1);

ntol=1;
for i=1:dim
    tmp{:,i}=Key1dMesh(Id{i},:);
    n(i)=size(tmp{:,i},1);
    ntol=ntol*n(i);
end
ntol
% key=zeros(ntol,3*dim);
count=1;

for dd=dim:-1:1
    for i=1:n(dd)
        key(count,[dd:dim:3*dim])=tmp{:,dd}(i,:);
        count=count+1;
    end
end

% key(:,3*dim-dim+1:3*dim)=key(:,3*dim-dim+1:3*dim)-1;

end
%
% function [key,count]=recurseloop(Id,count)
% dim=size(Id,1);
% if count=
%     key(count,[dd])
%     count=count+1;
% end
