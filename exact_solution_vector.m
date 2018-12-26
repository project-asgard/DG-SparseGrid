function fval = exact_solution_vector(HASHInv,pde,time)

% Returns the wavelet transformed exact solution

LevX = pde.dimensions{1}.lev;
LevV = pde.dimensions{2}.lev;
Deg  = pde.dimensions{1}.deg;

fx1 = forwardMWT(LevX,Deg,pde.params.Lmin,pde.params.Lmax,pde.ExactFx,pde.params);
fv1 = forwardMWT(LevV,Deg,pde.params.Vmin,pde.params.Vmax,pde.ExactFv,pde.params);
ft1 = pde.ExactFt(time);

fxList = {fx1};
fvList = {fv1};
ftList = {ft1};

fval = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);

end