function [fval] = combine_dimensions_2(fList,ftList,HASHInv,pde)

% NOTE: This is a 2D routine. For higher D we would need to generalise.
%
% Combine (via kron product) a set of 1D multiwavelet transforms to form
% the higher (here 2D) sparse-grid multiwavelet representation.

% fxList and fyList are a list of terms to be kron, then added together,
% i.e.,
%
% result = kron(fx1,fy1)*ft1 + kron(fx2,fy2)*ft2 + ...

Dim = pde.Dim;
Deg = pde.Deg;

N = numel(fList);

nHash = numel(HASHInv);

fval = sparse(Deg^Dim * nHash,1);

for j=1:N % loop over number of addative terms
    
    f = fList{j};
    
    if Dim==2
        
        f1 = f(:,1);
        f2 = f(:,2);
        
        ft = ftList{j};
        
        for i=1:nHash
            
            ll=HASHInv{i};
            
            I1=ll(3,1);
            I2=ll(3,2);
            
            
            index_I1=[(I1-1)*Deg+1 : I1*Deg];
            index_I2=[(I2-1)*Deg+1 : I2*Deg];
            
            Index = Deg^2*(i-1)+1:Deg^2*i;
            
            tmp = f2(index_I2)*ft*transpose(f1(index_I1));
            
            fval(Index,1) = fval(Index,1) + tmp(:);
            
        end
        
    end
    
end

end
