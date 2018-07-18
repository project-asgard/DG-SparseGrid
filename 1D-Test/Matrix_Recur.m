function A = Matrix_Recur(m,n,Deg)

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre( 1,Deg);

quad_num = 10;
[quad_x,quad_w]=lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);


load(['two_scale_rel_',num2str(Deg),'.mat'])

H0(find(abs(H0)<1e-5))=0;
G0(find(abs(G0)<1e-5))=0;

H1 = zeros(Deg);
G1 = zeros(Deg);

for j_x = 1:Deg
    for j_y = 1:Deg
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)  )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(Deg+j_x+j_y-2))*G0(j_x,j_y);
    end
end

H = [H0,H1];
G = [G0,G1];

% A11 = [Dp_val'*(quad_w.*p_val)]*2^0;
% A12 = H*kron(eye(2),A11)*G'*2;
% A21 = G*kron(eye(2),A11)*H'*2;
% A22 = G*kron(eye(2),A11)*G'*2;

A11 = -p_1'*p_1/2+p_2'*p_2/2+[Dp_val'*(quad_w.*p_val)];
A21 = -p_1'*p_2/2;
A12 =  p_2'*p_1/2;
% GradXFluxC = -blktridiag([Amd],[Asub],[Asup],nx);
A22 = G*kron(eye(2),A11)*G'*2;
% Adding Periodic Boundary Conditions
% % IndexEnd = dof_1D_x-k+1:dof_1D_x;
% % IndexSta = 1:k;
% % Iu = repmat(IndexEnd,k,1);
% % Iv = repmat(IndexSta,k,1);
% % GradXFluxC = GradXFluxC...
% %     +sparse([Iv',Iu'],[Iu,Iv],-[Asub,Asup],dof_1D_x,dof_1D_x);

if m==0 && n==0
    A = A11;
elseif m==0 && n==1
    A = A12;
elseif m==1 && n==0
    A = A21;
elseif m==1 && n==1
    A = A22;
else
    if m == 0
        tmp = Matrix_Recur(0,n-1,Deg);
        A = H*kron(tmp,eye(2))*2;
%         111
    elseif m==1
%         222
        tmp = Matrix_Recur(0,n-1,Deg);
        A = G*kron(tmp,eye(2))*2;
    elseif n==0
        tmp = Matrix_Recur(m-1,0,Deg);
        A = kron(eye(2),tmp)*H'*2;
    elseif n==1
        tmp = Matrix_Recur(m-1,0,Deg);
        A = kron(eye(2),tmp)*G'*2;
    else
        tmp = Matrix_Recur(m-1,n-1,Deg);
        A = kron(eye(2),tmp)*2;
    end
end

A = sparse(A);
end

