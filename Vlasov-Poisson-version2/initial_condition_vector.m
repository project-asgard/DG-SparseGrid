function fval = initial_condition_vector(fx,fv,Deg,Dim,HASHInv,pde)

nHash = numel(HASHInv);
Dim = pde.params.Dim;
Deg = pde.params.Deg;
DOF = Deg^Dim*nHash;

fxList = {fx};
fvList = {fv};
ftList = {zeros(DOF,1)+1};

fval = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);

end