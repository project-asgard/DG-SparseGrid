function v=dlegendre2(x,k,q)
% v=dlegendre2(x,k)
%
% q's Derivatives of Legendre Polynomials with degree k on [-1,1]
%
% P(x,0,0) = 1, P(x,1,0) = x
% P(x,0,1) = 0, P(x,1,1) = 1
%
% (n+1) * P(x,n+1,0) = (2*n+1)*x*P(x,n,0) - n * P(x,n-1,0)
%
% (n+1) * P(x,n+1,q) = q*(2*n+1)*P(x,n,q-1)+(2*n+1)*x*P(x,n,q)-n*P(x,n-1,q)
%
% integral( P(x,m) * P(x,n), x in [-1,1]) = 2/(2*n+1)*delta(m,n)
%
% ------------------------------------------------
nx = prod(size(x));
x = reshape(x, nx,1);
v = zeros(nx,k,q+1);

% Index starts from 1
% P1 -- > L0
% P2 -- > L1
maxn = k-1;


v(1:nx,1,1) = 1; % P(x,0,0) = 1

if (k >= 2),
    v(1:nx,2,1) = x; % P(x,1,0) = x
    
end;
if (q >= 1)
    v(1:nx,1,2) = 0;
    v(1:nx,2,2) = 1; % P(x,1,0) = x
end

% compute the functions value
if (k >= 3),
    Pnm1 = v(1:nx,1);
    Pn = v(1:nx,2);
    for n=1:(maxn-1),
        np1 = n + 1;
        Pnp1(1:nx) = ((2*n+1)*x(1:nx).*Pn(1:nx) - n * Pnm1(1:nx) )/(n+1);
        v(1:nx,1+np1,1) = Pnp1(1:nx);
        Pnm1(1:nx) = Pn(1:nx);
        Pn(1:nx) = Pnp1(1:nx);
    end;
end;

% compute the derivatives
if (k >= 3 && q >= 1),
    %     Pnm1 = v(1:nx,1,1); % P0
    %     DPnm1 = v(1:nx,1,2); % P0'
    %     Pn = v(1:nx,2,1); % P1
    %     DPn = v(1:nx,2,2); % P1'
    
    for n=1:(maxn-1),
        np1 = n + 1; % P2
        
        for p = 2:min(q+1,np1+1) % take derivatives for the Polynomial Pn
            % update
            DPnm1(1:nx,1) = v(1:nx,np1-1,p);
            DPn(1:nx,1) = v(1:nx,np1,p);
            Pn(1:nx,1) = v(1:nx,np1,p-1);
            
            DPnp1(1:nx) = ((p-1)*(2*n+1)*Pn(1:nx)...
                +(2*n+1)*x(1:nx).*DPn(1:nx) - n * DPnm1(1:nx) )/(n+1);
            v(1:nx,1+np1,p) = DPnp1(1:nx);
            
            
            % (n+1) * P(x,n+1,q) = q*(2*n+1)*P(x,n,q-1)+(2*n+1)*x*Pn(x,n,q)-n*Pnm1(x,n-1,q)
            
            
        end
        
    end;
end;



% -------------
% compute  norm
% -------------
norm2 = zeros(k,1);
for n=0:maxn,
    norm2(1 + n) = 2/(2*n+1);
end;

% ------------------------------
% normalize legendre polynomials
% ------------------------------
for n=0:maxn,
    dscale = 1/sqrt( norm2(1+n) );
    v(1:nx, 1+n,:) = v(1:nx,1+n,:) * dscale;
end;

% ----------------------------
% zero out points out of range
% ----------------------------
out_of_range = find( (x < -1) | (x > 1) );
v( out_of_range, 1:k) = 0;


% ----------------------------------------
% scaling to use  normalization, <Pn(x), Pn(x)> = 2
% ----------------------------------------
v = v * sqrt(2);


