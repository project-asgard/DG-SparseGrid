clc

Lev = 5;
Dim = 2;
Deg = 1;

count = 0;
DoFs = Deg^2*(Start_LevCell(0,Lev)+Num4Cell(0,Lev));
Mat_GL = sparse(DoFs,DoFs);

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
            
%             [Ix,Iy,Jx,Jy]
        
            tmpMat1 = Matrix_Recur(Ix,Jx,Deg);
            tmpMat2 = Matrix_Recur(Iy,Jy,Deg);
            
            IndexI_loc = GlobalIndexMap(Ix,Iy,Deg);
            IndexJ_loc = GlobalIndexMap(Jx,Jy,Deg);
            
            IndexI = IndexI_loc'*ones(1,numel(IndexJ_loc));
            IndexJ = ones(numel(IndexI_loc),1)*IndexJ_loc;
            
            num_size = length(IndexI);
            Mat_GL = Mat_GL+sparse(IndexI,IndexJ,kron(tmpMat1,tmpMat2),DoFs,DoFs);
            
            count = count+1;
            
            A_encode{count}.IndexI = IndexI;
            A_encode{count}.IndexJ = IndexJ;
            A_encode{count}.A1=tmpMat1;
            A_encode{count}.A2=tmpMat2;
            
             
        end
        
        
        
    end
end
