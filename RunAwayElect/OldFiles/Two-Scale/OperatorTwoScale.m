function FMWT_COMP = OperatorTwoScale(k,Np)
%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: k
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************
load(['two_scale_rel_',num2str(k),'.mat'])
n=log2(Np);


H0(find(abs(H0)<1e-5))=0;
G0(find(abs(G0)<1e-5))=0;

H1 = zeros(k);
G1 = zeros(k);

for j_x = 1:k
    for j_y = 1:k
        H1(j_x,j_y) = ((-1)^(j_x+j_y-2)  )*H0(j_x,j_y);
        G1(j_x,j_y) = ((-1)^(k+j_x+j_y-2))*G0(j_x,j_y);
    end
end

FMWT = zeros(k*Np);
iFMWT = zeros(k*Np);

for j=1:Np/2
    % The reverse order from Lin
    FMWT(k*(j-1)+1:k*j,2*k*(j-1)+1:2*k*j)=[H0 H1];
    FMWT(k*(j+Np/2-1)+1:k*(j+Np/2),2*k*(j-1)+1:2*k*j) = [G0 G1];
end
iFMWT=FMWT';

sp = [];
FMWT_COMP = eye(k*Np);
for j=1:n
    cFMWT = FMWT;
    % Reverse the index in matrix from Lin
    if j>1
        cFMWT = zeros(k*Np);
        cn = 2^(n-j+1)*k;
        cnr=Np*k-cn;
        cFMWT(cn+1:k*Np,cn+1:k*Np)=eye(Np*k-cn);
        cFMWT(1:cn/2,1:cn)=FMWT(1:cn/2,1:cn);
        cFMWT(cn/2+1:cn,1:cn)=FMWT(k*Np/2+1:k*Np/2+cn/2,1:cn);
    end

    FMWT_COMP = cFMWT*FMWT_COMP;
end

iFMWT_COMP=FMWT_COMP';
