function fval = initial_condition_vector(fx,fv,Deg,Dim,HASHInv,pde)

fxList = {fx};
fvList = {fv};
ftList = {1};

fval = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);

end