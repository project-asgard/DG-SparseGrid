function fval = initial_condition_vector(fx,fv,Deg,Dim,HASHInv)

nHash = numel(HASHInv);

fval = sparse(Deg^Dim * nHash,1);

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
            vv = kron( fv(Index1), fx(Index2) );
            mm = Deg^Dim * nHash;
            nn = 1;
            
            fval = fval + sparse(ii,jj,vv,mm,nn);
            
        end
    end  
end

end