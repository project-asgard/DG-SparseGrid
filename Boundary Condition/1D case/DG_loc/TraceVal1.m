function [avg,jump]=TraceVal1(x,k)
%------------------------------------------------
% compute the value on (-1)
% compute the value on the middle point (0)
% compute the value on (1)
%------------------------------------------------
jump=zeros(3,k);
avg=zeros(3,k);

f1=legendre(x,k);
d1=dlegendre(x,k);

ix1=find(abs(x+1)<1e-5);
ix2=find(abs(x-1)<1e-5);
ix=[ix1;ix2];

for i=1:k
% %     jump(ix1,i) = f1(ix1,i);
% %     jump(ix2,i) =-f1(ix2,i);
% %     
% %     avg(:,i)  = d1(:,i);
% %     avg(ix1,i)= d1(ix1,i)/2;
% %     avg(ix2,i)= d1(ix2,i)/2;

    jump(ix1,i) =-f1(ix1,i); % point -1
    jump(ix2,i) = f1(ix2,i); % point 1
    
    avg(:,i)  = d1(:,i);
    avg(ix1,i)= d1(ix1,i)/2;
    avg(ix2,i)= d1(ix2,i)/2;
end

avg=avg*2;


end

function v=legendre(x,k)
% Legendre Polynomials with degree k on [-1,1]
v(:,1)=x-x+1;
v(:,2)=x*sqrt(2*2-1);
v(:,3)=1/2*(3*x.^2-1)*sqrt(3*2-1);
v(:,4)=1/2*(5*x.^3-3*x)*sqrt(4*2-1);
v(:,5)=1/8*(35*x.^4-30*x.^2+3)*sqrt(5*2-1);
v(:,6)=1/8*(63*x.^5-70*x.^3+15*x)*sqrt(6*2-1);
v(:,7)=1/16*(231*x.^6-315*x.^4+105*x.^2-5)*sqrt(7*2-1);

v=v(:,1:k);
end

function v=dlegendre(x,k)
% Legendre Polynomials with degree k on [-1,1]
v(:,1)=x-x;
v(:,2)=x-x+1*sqrt(2*2-1);
v(:,3)=1/2*(3*2*x)*sqrt(3*2-1);
v(:,4)=1/2*(5*3*x.^2-3)*sqrt(4*2-1);
v(:,5)=1/8*(35*4*x.^3-30*2*x)*sqrt(5*2-1);
v(:,6)=1/8*(63*5*x.^4-70*3*x.^2+15)*sqrt(6*2-1);
v(:,7)=1/16*(231*6*x.^5-315*4*x.^3+105*2*x)*sqrt(7*2-1);

v=v(:,1:k);
end