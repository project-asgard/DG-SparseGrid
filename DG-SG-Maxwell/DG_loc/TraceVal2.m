function [avg,jump]=TraceVal2(x,k)
%--------------------------------
% compute the value on (-1)
% compute the value on the middle point (0)
% compute the value on (1)
%--------------------------------
x=2*x-1;

jump=zeros(3,k);
avg=zeros(3,k);

f1=f_basis(x,k);
d1=df_basis(x,k);


ix1=find(abs(x+1)<1e-5);
ix0=find(abs(x)<1e-5);
ix2=find(abs(x-1)<1e-5);


ix=find(x<0);
f1(ix,:)=f_basis(-x(ix),k);
d1(ix,:)=df_basis(-x(ix),k);


for i=1:k
%     jump(ix1,i) = (-1)^(i+k-1)*f1(ix1,i); % point -1
%     jump(ix2,i)= -f1(ix2,i); % point 1
%     jump(ix0,i)= (-(-1)^(i+k-1)+1)*f1(ix0,i); % point 0
%     
%     avg(:,i) = d1(:,i);
%     avg(ix,i)= -(-1)^(i+k-1)*d1(ix,i);
%     
%     avg(ix1,i)  = -(-1)^(i+k-1)*d1(ix1,i)/2;
%     avg(ix0,i)  = (1-(-1)^(i+k-1))*d1(ix0,i)/2;
%     avg(ix2,i)  = d1(ix2,i)/2;

    jump(ix1,i) =- (-1)^(i+k-1)*f1(ix1,i); % point -1
    jump(ix2,i)=  f1(ix2,i); % point 1
    jump(ix0,i)= ((-1)^(i+k-1)-1)*f1(ix0,i); % point 0
    
    avg(:,i) = d1(:,i);
    avg(ix,i)= -(-1)^(i+k-1)*d1(ix,i);
    
    avg(ix1,i)  = -(-1)^(i+k-1)*d1(ix1,i)/2;
    avg(ix0,i)  =  (1-(-1)^(i+k-1))*d1(ix0,i)/2;
    avg(ix2,i)  = d1(ix2,i)/2;
end


avg=2*avg;

end

function f=f_basis(x,k)
switch k
    case 1
        f(:,1)=x-x+sqrt(1/2);
    case 2
%         disp('case 2')
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