clc

Lev = 3;
Dim = 2;
Deg = 2;

count = 0;

for i_val = 0:Lev
    II = perm_eq(Dim,i_val);
    n1 = size(II,1);
   
    for j_val = 0:Lev
        JJ = perm_eq(Dim,j_val);
        n2 = size(JJ,1);
        JJ = repmat(JJ,n1,1);
                
        for i = 1:n1*n2
            Ix = II(ceil(i/n2),1);
            Iy = II(ceil(i/n2),2);
            Jx = JJ(i,1);
            Jy = JJ(i,2);
             
            [Ix,Iy,Jx,Jy]
        
            Mat1 = Matrix_Recur(Ix,Jx,Deg);
            Mat2 = Matrix_Recur(Iy,Jy,Deg);
            
            Index_I = GlobalIndexMap(Ix,Iy,Deg); 
            Index_J = GlobalIndexMap(Jx,Jy,Deg); 
            
            count = count+1; 
        end
        
        
        
    end
end

count