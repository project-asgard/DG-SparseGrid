function [fval] = combine_dimensions_2(fxList,fvList,ftList,HASHInv,pde)

% NOTE: This is a 2D routine. For higher D we would need to generalise.
%
% Combine (via kron product) a set of 1D multiwavelet transforms to form
% the higher (here 2D) sparse-grid multiwavelet representation.

% fxList and fyList are a list of terms to be kron, then added together,
% i.e.,
%
% result = kron(fx1,fy1)*ft1 + kron(fx2,fy2)*ft2 + ...

Dim = 2;

Deg = pde.params.Deg;

N = numel(fxList);

nHash = numel(HASHInv);

fval = sparse(Deg^Dim * nHash,1);



for j=1:N % loop over number of addative terms
    
    fx = fxList{j};
    fv = fvList{j};
    ft = ftList{j};
    
    for i=1:nHash
        
        ll=HASHInv{i};
        
        % 1D indices for (Lev1,Cell1)-->Index1,(Lev2,Cell2)-->Index2
        I1=ll(5);
        I2=ll(6);
        
        for k1 = 1:Deg
            
            Index1 = Deg*(I1-1)+k1;
            
            for k2 = 1:Deg
                
                Index2 = Deg*(I2-1)+k2;
                Index0 = Deg^2*(i-1)+Deg*(k1-1)+k2;
                
                ii = Index0;
                jj = 1;
                vv = kron( fv(Index1), fx(Index2) ) * ft;

                mm = Deg^Dim * nHash;
                nn = 1;
                
                fval = fval + sparse(ii,jj,vv,mm,nn);
                
            end
        end
    end
    
end

end