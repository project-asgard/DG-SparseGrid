function h=waveletbasis(x,k)
% basis defined on (0,1)

h = 2^(1/2)*orthbasis(2*x-1,k);

ix=find(x>1);
h(ix,:)=0;
ix=find(x<0);
h(ix,:)=0;

end

function f = orthbasis(x,k)
% wavelet functions f_1,\dots, f_k
% Ref:B. Alpert, A class of bases in L^2 for the sparse representation of
% integral operators, SIAM J. MATH. ANAL, 24(1993): 246-262. (Table 1)
% function f is defined in the interval (0,1) and then extended to the
% interval (-1,1), as f_i(x)=(-1)^{i+k-1}f_i(-x)
% first get the value for x\in (0,1)
ix=find(x>=0);
f(ix,:)=f_basis(x(ix),k);
ix=find(x<0);
tmp=f_basis(-x(ix),k);

ix0=find(abs(x)<1e-7);
ix1=find(abs(x-1)<1e-7);
ix2=find(abs(x+1)<1e-7);

f0=f_basis(0,k);
for i=1:k
    f(ix,i)=(-1)^(i+k-1)*tmp(:,i);
    
    if size(ix0,1)>0
        f(ix0(1),i)=(-1)^(i+k-1)*f0(i);
        f(ix0(2),i)=f0(i);
    end
    
end

if size(ix1,1)>1
    f(ix1(2),:)=0;
end

if size(ix2,1)>1
    f(ix2(1),:)=0;
end

end

function f=f_basis(x,k)
switch k
    case 1
        f(:,1)=x-x+sqrt(1/2);
    case 2
        f(:,1)=sqrt(3/2)*(-1+2*x);
        f(:,2)=sqrt(1/2)*(-2+3*x);
    case 3
        f(:,1)=1/3*sqrt(1/2)*(1-24*x+30*x.^2);
        f(:,2)=1/2*sqrt(3/2)*(3-16*x+15*x.^2);
        f(:,3)=1/3*sqrt(5/2)*(4-15*x+12*x.^2);
    case 4
        f(:,1)=sqrt(15/34)*(1+4*x-30*x.^2+28*x.^3);
        f(:,2)=sqrt(1/42)*(-4+105*x-300*x.^2+210*x.^3);
        f(:,3)=1/2*sqrt(35/34)*(-5+48*x-105*x.^2+64*x.^3);
        f(:,4)=1/2*sqrt(5/42)*(-16+105*x-192*x.^2+105*x.^3);
    case 5
        f(:,1)=sqrt(1/186)*(1+30*x+210*x.^2-840*x.^3+630*x.^4);
        f(:,2)=1/2*sqrt(1/38)*(-5-144*x+1155*x.^2-2240*x.^3+1260*x.^4);
        f(:,3)=sqrt(35/14694)*(22-735*x+3504*x.^2-5460*x.^3+2700*x.^4);
        f(:,4)=1/8*sqrt(21/38)*(35-512*x+1890*x.^2-2560*x.^3+1155*x.^4);
        f(:,5)=1/2*sqrt(7/158)*(32-315*x+960*x.^2-1155*x.^3+480*x.^4);
end

f=f(:,1:k);



end
