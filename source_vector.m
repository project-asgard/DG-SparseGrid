function fval = source_vector(HASHInv,pde,time)

% Returns the wavelet transformed source

LevX = pde.Lev(1);%params.LevX;
LevV = pde.Lev(2);%params.LevV;
Deg = pde.params.Deg;


Lmin = pde.params.Lmin;
Lmax = pde.params.Lmax;
Vmin = pde.params.Vmin;
Vmax = pde.params.Vmax;

fx1 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source1x,pde.params);
fv1 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source1v,pde.params);
ft1 = pde.source1t(time);

fx2 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source2x,pde.params);
fv2 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source2v,pde.params);
ft2 = pde.source2t(time);

fx3 = forwardMWT(LevX,Deg,Lmin,Lmax,pde.source3x,pde.params);
fv3 = forwardMWT(LevV,Deg,Vmin,Vmax,pde.source3v,pde.params);
ft3 = pde.source3t(time);

fxList = {fx1,fx2,fx3};
fvList = {fv1,fv2,fv3};
ftList = {ft1,ft2,ft3};

fval = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);

end