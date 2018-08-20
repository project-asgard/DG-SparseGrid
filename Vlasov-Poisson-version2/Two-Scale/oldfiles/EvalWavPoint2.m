function [f] = EvalWavPoint2(Lstart,Lend,maxLev,Deg,fcoef,x)
% This code evaluates the wavelet functions at given points
%
% Meval denotes the level and cel, where contains point x
% f denotes the value of func_wavelet(x)
%----------------------------------------------------------------------
load(['two_scale_rel_',num2str(Deg),'.mat']);

Lmax = Lend-Lstart;
Meval = zeros(maxLev+1,2);
MIndex = zeros(maxLev+1,1);
MIndex_full = zeros (Deg*(maxLev+1),1);

f = zeros(Deg*(maxLev+1),1);

% Lev = 0
hx = Lmax;
MidPoint = (Lend-Lstart)/2;
xhat = (x-MidPoint)*2/hx;
MIndex(1) = 1;
for k = 1:Deg
    val = polyval(scale_co(k,:),xhat);
    f(k) = val/sqrt(hx)*sqrt(Lmax/2);
    
    MIndex_full(k) = k;
end
Meval(1,1)=1;
% Lev = 1 : maxLev
for L = 1:maxLev
%     hx = Lmax/2^L;
    hx = Lmax/2^(L-1);
    if abs(x-Lstart)<1e-5
        cel = 1;
%     elseif abs(x-Lend)<1e-5
%         cel = 2^(L-1);
    else
        cel = ceil((x-Lstart)/hx);
    end
     MidPoint = ( (Lstart+cel*hx)+(Lstart+(cel-1)*hx) )/2;
     
    Meval(L+1,1)= cel;
    MIndex(L+1) = 2^(L-1)+cel;
    MIndex_full(Deg*L+[1:Deg]) = (2^(L-1))*Deg+Deg*(cel-1)+[1:Deg];
    % This turned on cell is from [hx*cel,hx*(cel+1)]
    % the middle point is hx*cel+hx/2
    % mapping the point from physical domain to reference domain
    xhat = (x-MidPoint)*2/hx;
    for k = 1:Deg
        if xhat<=0
            val = polyval(phi_co(k,:),xhat);
        else
            val = polyval(phi_co(k+Deg,:),xhat);
        end
        f(L*Deg+k) = val*1/sqrt(hx)*1/sqrt(Lmax);%*sqrt(2*k-1);
    end
    
end

% [Deg*(MIndex-1)+1 Deg*MIndex]
% 
% MIndex_full

f_loc = (sparse(MIndex_full,ones(Deg*(maxLev+1),1),f,Deg*2^(maxLev),1));

f = fcoef'*f_loc;
end