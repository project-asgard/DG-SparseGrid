function h=dwaveletbasis(x,k)
% basis defined on (0,1)
h = 2*2^(1/2)*dorthbasis (2*x-1,k);
end

function f = dorthbasis(x,k)
% wavelet functions f_1,\dots, f_k
% Ref:B. Alpert, A class of bases in L^2 for the sparse representation of
% integral operators, SIAM J. MATH. ANAL, 24(1993): 246-262. (Table 1)
% function f is defined in the interval (0,1) and then extended to the
% interval (-1,1), as f_i(x)=(-1)^{i+k-1}f_i(-x)
% first get the value for x\in (0,1)
f=zeros(size(x,1),k);
ix=find(x>=0);

f(ix,:)=df_basis(x(ix),k);
ix=find(x<0);
tmp=df_basis(-x(ix),k);
for i=1:k
f(ix,i)=(-1)^(i+k)*tmp(:,i);
end
end

function f=df_basis(x,k)

switch k
    case 1
        f(:,1)=x-x;
    case 2
        f(:,1)=sqrt(3/2)*(2)+x-x;
        f(:,2)=sqrt(1/2)*(3)+x-x;
    case 3
        f(:,1)=1/3*sqrt(1/2)*(-24+30*2*x);
        f(:,2)=1/2*sqrt(3/2)*(-16+15*2*x);
        f(:,3)=1/3*sqrt(5/2)*(-15+12*2*x);
    case 4
        f(:,1)=sqrt(15/34)*(4-30*2*x+28*3*x.^2);
        f(:,2)=sqrt(1/42)*(105-300*2*x+210*3*x.^2);
        f(:,3)=1/2*sqrt(35/34)*(48-105*2*x+64*3*x.^2);
        f(:,4)=1/2*sqrt(5/42)*(105-192*2*x+105*3*x.^2);
    case 5
        f(:,1)=sqrt(1/186)*(30+210*2*x-840*3*x.^2+630*4*x.^3);
        f(:,2)=1/2*sqrt(1/38)*(-144+1155*2*x-2240*3*x.^2+1260*4*x.^3);
        f(:,3)=sqrt(35/14694)*(-735+3504*2*x-5460*3*x.^2+2700*4*x.^3);
        f(:,4)=1/8*sqrt(21/38)*(-512+1890*2*x-2560*3*x.^2+1155*4*x.^3);
        f(:,5)=1/2*sqrt(7/158)*(-315+960*2*x-1155*3*x.^2+480*4*x.^3);
end
f=f(:,1:k);

end
