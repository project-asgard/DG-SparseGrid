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

                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                
                A_encode{count}.B1=M_mass(I1,I1);
                A_encode{count}.B2=Stiff_1D(I2,I2);
                A_encode{count}.B3=M_mass(I3,I3);
                A_encode{count}.B4=M_mass(I4,I4);

                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.B1, ...
                                  A_encode{count}.B2, ...
                                  A_encode{count}.B3, ...
                                  A_encode{count}.B4 );
                
                A_encode{count}.C1=M_mass(I1,I1);
                A_encode{count}.C2=M_mass(I2,I2);
                A_encode{count}.C3=Stiff_1D(I3,I3);
                A_encode{count}.C4=M_mass(I4,I4);

                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.C1, ...
                                  A_encode{count}.C2, ...
                                  A_encode{count}.C3, ...
                                  A_encode{count}.C4 );
                
                A_encode{count}.D1=M_mass(I1,I1);
                A_encode{count}.D2=M_mass(I2,I2);
                A_encode{count}.D3=M_mass(I3,I3);
                A_encode{count}.D4=Stiff_1D(I4,I4);

                kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.D1, ...
                                  A_encode{count}.D2, ...
                                  A_encode{count}.D3, ...
                                  A_encode{count}.D4 );
                
                A_encode{count}.IndexI=double(Index_I);
                A_encode{count}.IndexJ=double(Index_I);

                kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                
                count=count+1;
                
                % Term 2: S*I*I*I--Assume j2_level==i2_level
                % j3_level==i3_level j4_level==i4_level
                j2_level=i2_level;j3_level=i3_level;j4_level=i4_level;
                for j1_level=0:i1_level-1
                    
                    J1=Index_1D(k,j1_level);
                    J2=I2;J3=I3;J4=I4;
                    
                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);
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

                    
                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.A1, ...
                                  A_encode{count}.A2, ...
                                  A_encode{count}.A3, ...
                                  A_encode{count}.A4 );
                    
                    A_encode{count}.B1=0;
                    A_encode{count}.B2=0;
                    A_encode{count}.B3=0;
                    A_encode{count}.B4=0;
                    
                    A_encode{count}.C1=0;
                    A_encode{count}.C2=0;
                    A_encode{count}.C3=0;
                    A_encode{count}.C4=0;
                    
                    A_encode{count}.D1=0;
                    A_encode{count}.D2=0;
                    A_encode{count}.D3=0;
                    A_encode{count}.D4=0;
                    
                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);

                    
                    A_encode{count+1}.A1=Stiff_1D(J1,I1)';
                    A_encode{count+1}.A2=M_mass(J2,I2)';
                    A_encode{count+1}.A3=M_mass(J3,I3)';
                    A_encode{count+1}.A4=M_mass(J4,I4)';

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count+1}.A1, ...
                                  A_encode{count+1}.A2, ...
                                  A_encode{count+1}.A3, ...
                                  A_encode{count+1}.A4 );
                    
                    A_encode{count+1}.B1=0;
                    A_encode{count+1}.B2=0;
                    A_encode{count+1}.B3=0;
                    A_encode{count+1}.B4=0;
                    
                    A_encode{count+1}.C1=0;
                    A_encode{count+1}.C2=0;
                    A_encode{count+1}.C3=0;
                    A_encode{count+1}.C4=0;
                    
                    A_encode{count+1}.D1=0;
                    A_encode{count+1}.D2=0;
                    A_encode{count+1}.D3=0;
                    A_encode{count+1}.D4=0;
                    
                    A_encode{count+1}.IndexI=Index_I;
                    A_encode{count+1}.IndexJ=Index_J;

                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexJ);
                    
                    count=count+2;
                    clear Index_J
                end
                
                % Term 3: I*S*I*I--Assume j1_level==i1_level
                % j3_level==i3_level j4_level==i4_level
                j1_level=i1_level;j3_level=i3_level;j4_level=i4_level;
                for j2_level=0:i2_level-1
                    
                    J2=Index_1D(k,j2_level);
                    J1=I1;J3=I3;J4=I4;
                    
                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);
                    for jjj=1:size(key_j,1)
                        Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
                    end
                    
                    
                    
                    [II,JJ]=meshgrid(Index_J,Index_I);
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),Stiff_1D(J2,I2)),M_mass(J3,I3)),M_mass(J4,I4));
                    
                    A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=0;
                    A_encode{count}.A2=0;
                    A_encode{count}.A3=0;
                    A_encode{count}.A4=0;
                    
                    A_encode{count}.B1=M_mass(J1,I1);
                    A_encode{count}.B2=Stiff_1D(J2,I2);
                    A_encode{count}.B3=M_mass(J3,I3);
                    A_encode{count}.B4=M_mass(J4,I4);

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.B1, ...
                                  A_encode{count}.B2, ...
                                  A_encode{count}.B3, ...
                                  A_encode{count}.B4 );
                    
                    A_encode{count}.C1=0;
                    A_encode{count}.C2=0;
                    A_encode{count}.C3=0;
                    A_encode{count}.C4=0;
                    
                    A_encode{count}.D1=0;
                    A_encode{count}.D2=0;
                    A_encode{count}.D3=0;
                    A_encode{count}.D4=0;
                    
                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    
                    A_encode{count+1}.A1=0;
                    A_encode{count+1}.A2=0;
                    A_encode{count+1}.A3=0;
                    A_encode{count+1}.A4=0;
                    
                    A_encode{count+1}.B1=M_mass(J1,I1)';
                    A_encode{count+1}.B2=Stiff_1D(J2,I2)';
                    A_encode{count+1}.B3=M_mass(J3,I3)';
                    A_encode{count+1}.B4=M_mass(J4,I4)';

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count+1}.B1, ...
                                  A_encode{count+1}.B2, ...
                                  A_encode{count+1}.B3, ...
                                  A_encode{count+1}.B4 );
                    
                    A_encode{count+1}.C1=0;
                    A_encode{count+1}.C2=0;
                    A_encode{count+1}.C3=0;
                    A_encode{count+1}.C4=0;
                    
                    A_encode{count+1}.D1=0;
                    A_encode{count+1}.D2=0;
                    A_encode{count+1}.D3=0;
                    A_encode{count+1}.D4=0;
                    
                    A_encode{count+1}.IndexI=Index_I;
                    A_encode{count+1}.IndexJ=Index_J;

                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexJ);

                    count=count+2;
                    clear Index_J
                end
                
                % Term 4: I*I*S*I--Assume j1_level==i1_level
                % j2_level==i2_level j4_level==i4_level
                j1_level=i1_level;j2_level=i2_level;j4_level=i4_level;
                for j3_level=0:i3_level-1
                    
                    J3=Index_1D(k,j3_level);
                    J1=I1;J2=I2;J4=I4;
                    
                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);
                    for jjj=1:size(key_j,1)
                        Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
                    end
                    
                    
                    
                    [II,JJ]=meshgrid(Index_J,Index_I);
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),M_mass(J2,I2)),Stiff_1D(J3,I3)),M_mass(J4,I4));
                    A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=0;
                    A_encode{count}.A2=0;
                    A_encode{count}.A3=0;
                    A_encode{count}.A4=0;
                    
                    A_encode{count}.B1=0;
                    A_encode{count}.B2=0;
                    A_encode{count}.B3=0;
                    A_encode{count}.B4=0;
                    
                    A_encode{count}.C1=M_mass(J1,I1);
                    A_encode{count}.C2=M_mass(J2,I2);
                    A_encode{count}.C3=Stiff_1D(J3,I3);
                    A_encode{count}.C4=M_mass(J4,I4);

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.C1, ...
                                  A_encode{count}.C2, ...
                                  A_encode{count}.C3, ...
                                  A_encode{count}.C4 );
                    
                    A_encode{count}.D1=0;
                    A_encode{count}.D2=0;
                    A_encode{count}.D3=0;
                    A_encode{count}.D4=0;
                    
                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    
                    A_encode{count+1}.A1=0;
                    A_encode{count+1}.A2=0;
                    A_encode{count+1}.A3=0;
                    A_encode{count+1}.A4=0;
                    
                    A_encode{count+1}.B1=0;
                    A_encode{count+1}.B2=0;
                    A_encode{count+1}.B3=0;
                    A_encode{count+1}.B4=0;
                    
                    A_encode{count+1}.C1=M_mass(J1,I1)';
                    A_encode{count+1}.C2=M_mass(J2,I2)';
                    A_encode{count+1}.C3=Stiff_1D(J3,I3)';
                    A_encode{count+1}.C4=M_mass(J4,I4)';

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count+1}.C1, ...
                                  A_encode{count+1}.C2, ...
                                  A_encode{count+1}.C3, ...
                                  A_encode{count+1}.C4 );
                    
                    A_encode{count+1}.D1=0;
                    A_encode{count+1}.D2=0;
                    A_encode{count+1}.D3=0;
                    A_encode{count+1}.D4=0;
                    
                    A_encode{count+1}.IndexI=Index_I;
                    A_encode{count+1}.IndexJ=Index_J;

                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexJ);
                    
                    count=count+2;
                    clear Index_J
                end
                
                % Term 5: I*I*I*S--Assume j1_level==i1_level
                % j2_level==i2_level j3_level==i3_level
                j1_level=i1_level;j2_level=i2_level;j3_level=i3_level;
                for j4_level=0:i4_level-1
                    
                    J4=Index_1D(k,j4_level);
                    J1=I1;J2=I2;J3=I3;
                    
                    key_j=GenerateKey4Dv2(J1,J2,J3,J4,Key1dMesh);
                    for jjj=1:size(key_j,1)
                        Index_J(jjj,1)=database.(sprintf('i%g_',key_j(jjj,:)));
                    end
 
                    [II,JJ]=meshgrid(Index_J,Index_I);
                    
                    tmp=kron(kron(kron(M_mass(J1,I1),M_mass(J2,I2)),M_mass(J3,I3)),Stiff_1D(J4,I4));
                    A_s=A_s+sparse([II,JJ],[JJ,II],[tmp',tmp'],dof_sparse,dof_sparse);
                    
                    % save matrices to files
                    A_encode{count}.A1=0;
                    A_encode{count}.A2=0;
                    A_encode{count}.A3=0;
                    A_encode{count}.A4=0;
                    
                    A_encode{count}.B1=0;
                    A_encode{count}.B2=0;
                    A_encode{count}.B3=0;
                    A_encode{count}.B4=0;
                    
                    A_encode{count}.C1=0;
                    A_encode{count}.C2=0;
                    A_encode{count}.C3=0;
                    A_encode{count}.C4=0;
                    
                    A_encode{count}.D1=M_mass(J1,I1);
                    A_encode{count}.D2=M_mass(J2,I2);
                    A_encode{count}.D3=M_mass(J3,I3);
                    A_encode{count}.D4=Stiff_1D(J4,I4);

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count}.D1, ...
                                  A_encode{count}.D2, ...
                                  A_encode{count}.D3, ...
                                  A_encode{count}.D4 );
                    
                    A_encode{count}.IndexI=Index_J;
                    A_encode{count}.IndexJ=Index_I;

                    kron_nnz = kron_nnz + length( A_encode{count}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count}.IndexJ);
                    
                    A_encode{count+1}.A1=0;
                    A_encode{count+1}.A2=0;
                    A_encode{count+1}.A3=0;
                    A_encode{count+1}.A4=0;
                    
                    A_encode{count+1}.B1=0;
                    A_encode{count+1}.B2=0;
                    A_encode{count+1}.B3=0;
                    A_encode{count+1}.B4=0;
                    
                    A_encode{count+1}.C1=0;
                    A_encode{count+1}.C2=0;
                    A_encode{count+1}.C3=0;
                    A_encode{count+1}.C4=0;
                    
                    A_encode{count+1}.D1=M_mass(J1,I1)';
                    A_encode{count+1}.D2=M_mass(J2,I2)';
                    A_encode{count+1}.D3=M_mass(J3,I3)';
                    A_encode{count+1}.D4=Stiff_1D(J4,I4)';

                    kron_flops = kron_flops + ...
                              kron_mult_cost4( ...
                                  A_encode{count+1}.D1, ...
                                  A_encode{count+1}.D2, ...
                                  A_encode{count+1}.D3, ...
                                  A_encode{count+1}.D4 );
                    
                    A_encode{count+1}.IndexI=Index_I;
                    A_encode{count+1}.IndexJ=Index_J;

                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexI);
                    kron_nnz = kron_nnz + length( A_encode{count+1}.IndexJ);
                    
                    count=count+2;
                    clear Index_J
                end
                
                clear Index_I
            end
        end
    end
    
end

save(['./Data/A_4D_encode.mat'],'A_encode');

% -----------------------
% ensure A_s is symmetric
% -----------------------
A_s  = (A_s + A_s')/2;

disp(sprintf('Np=%d,k=%d, n=%d, kron_flops=%g, kron_nnz=%g, nnz(A_s)=%g ', ...
              Np, k, size(A_s,1), ...
              kron_flops,    kron_nnz, nnz(A_s)  ));

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
  %celldisp (Aencode(1))
  N = prod(size(Aencode));
  y = zeros(size(x));
  for count = 1: N
    
    nkron = 4;
    IndexI = Aencode{count}.IndexI;
    IndexJ = Aencode{count}.IndexJ;
%    printf("Count = %d  size(indexI) = %d  %d  size(indexJ) = %d  %d\n", ...
%            count, size(IndexI), size(IndexJ));
    my_X_I = x(IndexI);
    my_X_J = x(IndexJ);
    
    do_transpose = not(all(size(IndexI) == size(IndexJ)) && all(IndexI == IndexJ))
    
%    printf("my_X_I = %f  \n", my_X_I);
%    printf("my_X_J = %f  \n", my_X_J);
    
    Acell{1} = Aencode{count}.A1';
    Acell{2} = Aencode{count}.A2';
    Acell{3} = Aencode{count}.A3';
    Acell{4} = Aencode{count}.A4';
    if (sum(sum(Acell{1})) && sum(sum(Acell{2})) &&  ...
        sum(sum(Acell{3})) && sum(sum(Acell{4}))),
      my_Y = kron_multd(nkron, Acell, my_X_I);
%      count
%      size(my_X)
%      size(my_Y)
      y(IndexJ) = y(IndexJ) + my_Y;

      if (do_transpose),
        Acell{1} = Aencode{count}.A1;
        Acell{2} = Aencode{count}.A2;
        Acell{3} = Aencode{count}.A3;
        Acell{4} = Aencode{count}.A4;
        my_Y = kron_multd(nkron, Acell, my_X_J);
        y(IndexI) = y(IndexI) + my_Y;
      endif;      
   endif;
    
    Bcell{1} = Aencode{count}.B1';
    Bcell{2} = Aencode{count}.B2';
    Bcell{3} = Aencode{count}.B3';
    Bcell{4} = Aencode{count}.B4';
    if (sum(sum(Bcell{1})) && sum(sum(Bcell{2})) &&  sum(sum(Bcell{3})) ...
      && sum(sum(Bcell{4}))),
      my_Y = kron_multd(nkron, Bcell, my_X_I);
      y(IndexJ) = y(IndexJ) + my_Y;
      
      if (do_transpose),
        Bcell{1} = Aencode{count}.B1;
        Bcell{2} = Aencode{count}.B2;
        Bcell{3} = Aencode{count}.B3;
        Bcell{4} = Aencode{count}.B4;
        my_Y = kron_multd(nkron, Bcell, my_X_J);
        y(IndexI) = y(IndexI) + my_Y;
      endif;
    endif;
    
    
    Ccell{1} = Aencode{count}.C1';
    Ccell{2} = Aencode{count}.C2';
    Ccell{3} = Aencode{count}.C3';
    Ccell{4} = Aencode{count}.C4';
    if (sum(sum(Ccell{1})) && sum(sum(Ccell{2})) &&  sum(sum(Ccell{3})) ...
      && sum(sum(Ccell{4}))),
      my_Y = kron_multd(nkron, Ccell, my_X_I);
      y(IndexJ) = y(IndexJ) + my_Y;

      if (do_transpose),
        Ccell{1} = Aencode{count}.C1;
        Ccell{2} = Aencode{count}.C2;
        Ccell{3} = Aencode{count}.C3;
        Ccell{4} = Aencode{count}.C4; 
        my_Y = kron_multd(nkron, Ccell, my_X_J);
        y(IndexI) = y(IndexI) + my_Y;
      endif;
    endif;    
    
    Dcell{1} = Aencode{count}.D1';
    Dcell{2} = Aencode{count}.D2';
    Dcell{3} = Aencode{count}.D3';
    Dcell{4} = Aencode{count}.D4';
    if (sum(sum(Dcell{1})) && sum(sum(Dcell{2})) &&  sum(sum(Dcell{3})) ...
      && sum(sum(Dcell{4}))),
      my_Y = kron_multd(nkron, Dcell, my_X_I);
      size(my_X_I);
      size(my_Y);
      size(y(IndexJ));
      y(IndexJ) = y(IndexJ) + my_Y;
      
      if (do_transpose),
        Dcell{1} = Aencode{count}.D1;
        Dcell{2} = Aencode{count}.D2;
        Dcell{3} = Aencode{count}.D3;
        Dcell{4} = Aencode{count}.D4;
        my_Y = kron_multd(nkron, Dcell, my_X_J);
        y(IndexI) = y(IndexI) + my_Y;
      endif;
    endif;    
  endfor 
endfunction

testX = ones(size(b_s));
Y1 = apply_As(testX, A_s);
Y2 = apply_Aencode(testX, A_encode);

size(Y1)
Y1(1:200)
size(Y2)
0.5 * Y2(1:200)

sum(Y1 - 0.5 * Y2)
exit

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
  use_Aencode = true
  time_iter = -time();
  if (use_function_handle),
    if (use_Aencode),
      myfun = "apply_Aencode";
      myarg = A_encode;
    else
      myfun = "apply_As"
      myarg = A_s
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
