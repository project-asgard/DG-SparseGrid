function ff = Rhs4CurlCurl(Deg,Lev,Lmax,pde)

Lend = Lmax;Lstart = 0;

quad_num = 10;
[quad_x,quad_w] = lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);


fx = zeros(quad_num,3);
fy = zeros(quad_num,3);
fz = zeros(quad_num,3);

nx = 2^Lev;
hx = (Lend-Lstart)/nx;


for i=0:2^Lev-1
    
    xi_x = hx*(quad_x/2+1/2+i)+Lstart;
    fx = pde.rhs(xi_x,1/2+xi_x-xi_x,1/2+xi_x-xi_x);
    fy = pde.rhs(1/2+xi_x-xi_x,xi_x,1/2+xi_x-xi_x);
    fz = pde.rhs(1/2+xi_x-xi_x,1/2+xi_x-xi_x,xi_x);
    
    
    index = Deg*i+1:Deg*(i+1);
    f1x(index,1) = p_val'*(quad_w.*fx(:,1))*hx*sqrt(1/hx)/2;
    f1y(index,1) = p_val'*(quad_w.*fy(:,1))*hx*sqrt(1/hx)/2;
    f1z(index,1) = p_val'*(quad_w.*fz(:,1))*hx*sqrt(1/hx)/2;
    %     f1 = kron(kron(f1x,f1y),f1z);
    
    
    f2x(index,1) = p_val'*(quad_w.*fx(:,2))*hx*sqrt(1/hx)/2;
    f2y(index,1) = p_val'*(quad_w.*fy(:,2))*hx*sqrt(1/hx)/2;
    f2z(index,1) = p_val'*(quad_w.*fz(:,2))*hx*sqrt(1/hx)/2;
    %     f2 = kron(kron(f2x,f2y),f2z);
    
    f3x(index,1) = p_val'*(quad_w.*fx(:,3))*hx*sqrt(1/hx)/2;
    f3y(index,1) = p_val'*(quad_w.*fy(:,3))*hx*sqrt(1/hx)/2;
    f3z(index,1) = p_val'*(quad_w.*fz(:,3))*hx*sqrt(1/hx)/2;
    %     f3 = kron(kron(f3x,f3y),f3z);
    
end

% convert to multiwavelet basis
FMWT_COMP = OperatorTwoScale(Deg,2^Lev);

f1x = FMWT_COMP*f1x;
f1y = FMWT_COMP*f1y;
f1z = FMWT_COMP*f1z;

f2x = FMWT_COMP*f2x;
f2y = FMWT_COMP*f2y;
f2z = FMWT_COMP*f2z;

f3x = FMWT_COMP*f3x;
f3y = FMWT_COMP*f3y;
f3z = FMWT_COMP*f3z;

ff = struct(...
    'f1x',f1x,'f1y',f1y,'f1z',f1z,...
    'f2x',f2x,'f2y',f2y,'f2z',f2z,...
    'f3x',f3x,'f3y',f3y,'f3z',f3z);