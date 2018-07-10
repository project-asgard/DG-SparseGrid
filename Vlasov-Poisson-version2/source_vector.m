function fval = source_vector(LevX,LevV,Deg,Lmax,Vmax,HASHInv,pde,time)
% projection of given source function to sparse grids
% here we assume source=source1+source2

fx1 = ProjCoef2Wav_v2(LevX,Deg,0, Lmax,pde.source1x);
fv1 = ProjCoef2Wav_v2(LevV,Deg,-Vmax,Vmax,pde.source1v);
ft1 = pde.source1t(time);

fx2 = ProjCoef2Wav_v2(LevX,Deg,0,Lmax,pde.source2x);
fv2 = ProjCoef2Wav_v2(LevV,Deg,-Vmax,Vmax,pde.source2v);
ft2 = pde.source2t(time);

nHash = numel(HASHInv);
Dim = 2;
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
            vv = kron( fv1(Index1), fx1(Index2) )*ft1+...
                 kron( fv2(Index1), fx2(Index2) )*ft2;
            mm = Deg^Dim * nHash;
            nn = 1;
            
            fval = fval + sparse(ii,jj,vv,mm,nn);
            
        end
    end  
end



end