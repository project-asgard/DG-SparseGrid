function v=lin_legendre(x,k)

% Legendre Polynomials with degree k on [-1,1]

use_matlab_legendre = true;
if use_matlab_legendre
    
    for kk=1:k
        n = kk-1;
        vk = legendre(n,x)*sqrt((n+1)*2-1);
        v(:,kk) = vk(1,:);
    end
    
else
    
    if k>7
        error('ERROR: max supported deg=9');
    end
    v(:,1)=x-x+1;
    v(:,2)=x*sqrt(2^2-1);
    v(:,3)=1/2*(3*x.^2-1)*sqrt(3*2-1);
    v(:,4)=1/2*(5*x.^3-3*x)*sqrt(4*2-1);
    v(:,5)=1/8*(35*x.^4-30*x.^2+3)*sqrt(5*2-1);
    v(:,6)=1/8*(63*x.^5-70*x.^3+15*x)*sqrt(6*2-1);
    v(:,7)=1/16*(231*x.^6-315*x.^4+105*x.^2-5)*sqrt(7*2-1);
    v=v(:,1:k);
    
end

ix=find(x>1);
v(ix,:)=0;
ix=find(x<-1);
v(ix,:)=0;