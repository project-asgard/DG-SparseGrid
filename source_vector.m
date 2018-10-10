function fval = source_vector(HASHInv,pde,time)

% Returns the wavelet transformed source

LevX = pde.Lev(1);%params.LevX;
LevV = pde.Lev(2);%params.LevV;
Deg = pde.Deg;
Dim = pde.Dim;
Lev = pde.Lev;
domain = pde.domain;

% Lmin = pde.params.Lmin;
% Lmax = pde.params.Lmax;
% Vmin = pde.params.Vmin;
% Vmax = pde.params.Vmax;
% 
% fx1 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source1x,pde.params);
% fv1 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source1v,pde.params);
ft1 = pde.source1t(time);

% fx2 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source2x,pde.params);
% fv2 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source2v,pde.params);
ft2 = pde.source2t(time);

% fx3 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source3x,pde.params);
% fv3 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source3v,pde.params);
ft3 = pde.source3t(time);

for d = 1:Dim
    f1(:,d) = forwardMWT(Lev(d),Deg,domain(1,d),domain(2,d),pde.source1{d},pde.params);
    f2(:,d) = forwardMWT(Lev(d),Deg,domain(1,d),domain(2,d),pde.source2{d},pde.params);
    f3(:,d) = forwardMWT(Lev(d),Deg,domain(1,d),domain(2,d),pde.source3{d},pde.params);
end

%fxList = {fx1,fx2,fx3};
%fvList = {fv1,fv2,fv3};
fList = {f1,f2,f3};
ftList = {ft1,ft2,ft3};

% fval = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);
fval = combine_dimensions_2(fList,ftList,HASHInv,pde);

end