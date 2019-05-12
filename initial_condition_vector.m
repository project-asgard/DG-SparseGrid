% function fval = initial_condition_vector(fx,fv,HASHInv,pde)
function fval = initial_condition_vector(HASHInv,pde,time)
%
% fxList = {fx};
% fvList = {fv};
% ftList = {1};

dims = pde.dimensions;
nDims = numel(dims);

for d=1:nDims
    fList{d} = forwardMWT(pde,dims{d}.lev,dims{d}.deg,dims{d}.domainMin,dims{d}.domainMax,dims{d}.init_cond_fn,pde.params);
end

% fx = pde.dimensions{1}.f0;
% fv = pde.dimensions{2}.f0;

% fval = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde);

ft = 1;
fval = combine_dimensions_D(fList,ft,HASHInv,pde);

end