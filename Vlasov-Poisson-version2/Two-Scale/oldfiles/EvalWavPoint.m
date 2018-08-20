function [f] = EvalWavPoint(Lstart,Lend,maxLev,Deg,fcoef,x)
% This code evaluates the wavelet functions at given points
%
% Meval denotes the level and cel, where contains point x
% f denotes the value of func_wavelet(x)
%----------------------------------------------------------------------
load(['two_scale_rel_',num2str(Deg),'.mat']);

Lmax = Lend-Lstart;
Meval = zeros(maxLev+1,2);
MIndex = zeros(maxLev+1,1);


nz = length(x);

MIndex_full = zeros (Deg*(maxLev+1),nz);
f = zeros(Deg*(maxLev+1),nz);

% Lev = 0
hx = Lmax;
MidPoint = (Lend-Lstart)/2;
xhat = (x-MidPoint)*2/hx;
MIndex(1) = 1;
for k = 1:Deg
    val = polyval(scale_co(k,:),xhat);
<<<<<<< HEAD
%     f(k) = val/sqrt(hx/2);%*1/sqrt(Lmax);
    f(k) = val/sqrt(hx)*1/sqrt(Lmax);
=======
    f(k,:) = val/sqrt(hx);%*sqrt(2*k-1);
>>>>>>> 22859c13fa5348f88c1989ae0cda2d5dca636de8
    
    MIndex_full(k,:) = k;
end

Meval(1,1)=1;
% Lev = 1 : maxLev
for L = 1:maxLev
%     hx = Lmax/2^L;
    hx = Lmax/2^(L-1);
    ix = find(abs(x-Lstart)<1e-5);
    cel = ceil((x-Lstart)/hx);
    cel(ix) = 1;

     MidPoint = ( (Lstart+cel*hx)+(Lstart+(cel-1)*hx) )/2;
     
% %     Meval(L+1,1)= cel;
% %     MIndex(L+1) = 2^(L-1)+cel;
    for k=1:Deg
    MIndex_full(Deg*L+k,1:nz) = (2^(L-1))*Deg+Deg*(cel'-1)+k;
    end
    % This turned on cell is from [hx*cel,hx*(cel+1)]
    % the middle point is hx*cel+hx/2
    % mapping the point from physical domain to reference domain
    xhat = (x-MidPoint)*2/hx;
    ip = find(xhat>0);
    
    for k = 1:Deg
%         if xhat<=0
            val = polyval(phi_co(k,:),xhat);
<<<<<<< HEAD
        else
            val = polyval(phi_co(k+Deg,:),xhat);
        end
%         hx = Lmax/2^L;
        f(L*Deg+k) = val/sqrt(hx)*1/sqrt(Lmax);
=======
%         else
            val(ip) = polyval(phi_co(k+Deg,:),xhat(ip));
%         end
        f(L*Deg+k,:) = val*1/sqrt(hx)*1/sqrt(Lmax);%*sqrt(2*k-1);
>>>>>>> 22859c13fa5348f88c1989ae0cda2d5dca636de8
    end
    
end

% [Deg*(MIndex-1)+1 Deg*MIndex]
% 
% MIndex_full

f_loc = (sparse(MIndex_full,ones(Deg*(maxLev+1),1)*[1:nz],f,Deg*2^(maxLev),nz));

f = f_loc'*fcoef;
end