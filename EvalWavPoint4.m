function [f_loc] = EvalWavPoint4(Lmin,Lmax,Lev,Deg,x)
% This code evaluates the wavelet functions at given points
%
% Meval denotes the level and cel, where contains point x
% f_loc denotes the value of func_wavelet(x)
%----------------------------------------------------------------------
load(['Two-Scale/two_scale_rel_',num2str(Deg),'.mat']);

range = Lmax-Lmin;

nz = length(x);

MIndex_full = zeros (Deg*(Lev+1),nz);
f = zeros(Deg*(Lev+1),nz);

% Lev = 0
hx = range;
MidPoint = (Lmax+Lmin)/2;
xhat = (x-MidPoint)*2/hx;
coef = 1/sqrt(hx);

for k = 1:Deg
    val = polyval(scale_co(k,:),xhat);
    f(k,:) = val*coef;
    
    MIndex_full(k,:) = k;
end

% Lev = 1 : Lev
for L = 1:Lev
    
    hx = range/2^(L-1);
    ix = find(abs(x-Lmin)<1e-5);
    cel = ceil((x-Lmin)/hx);
    cel(ix) = 1;
    
    MidPoint = ( (Lmin+cel*hx)+(Lmin+(cel-1)*hx) )/2;
    
    for k=1:Deg
        MIndex_full(Deg*L+k,1:nz) = (2^(L-1))*Deg+Deg*(cel'-1)+k;
    end
    % This turned on cell is from [hx*cel,hx*(cel+1)]
    % the middle point is hx*cel+hx/2
    % mapping the point from physical domain to reference domain
    xhat = (x-MidPoint)*2/hx;
    ip = find(xhat>0);
    
    for k = 1:Deg
        
        val = polyval(phi_co(k,:),xhat);
        
        val(ip) = polyval(phi_co(k+Deg,:),xhat(ip));
        
        coef = sqrt(1/hx);
        
        f(L*Deg+k,:) = val*coef;
    end
    
end

% f_loc evaluates DoF_1D basis at nz points xx
f_loc = (sparse(MIndex_full,ones(Deg*(Lev+1),1)*[1:nz],f,Deg*2^(Lev),nz));

end
