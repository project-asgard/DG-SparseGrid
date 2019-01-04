function Mat_tmp = GlobalMatrix(Mat,IHash,Deg)

DoFs = Deg*size(IHash,2);
Mat_tmp = sparse(DoFs,DoFs);

    for i = 1:size(IHash,2)
        for j = 1:size(IHash,2)
        l1 = IHash{i};
        n1 = l1(1);
        c1 = l1(2);
        loc = l1(3);%2^(n1-1)+c1+1;
        id1 = Deg*(loc-1)+[1:Deg];
        i1 = Deg*(i-1)+[1:Deg];
        
        
        l2 = IHash{j};
        n2 = l2(1);
        c2 = l2(2);
        loc = l2(3);%2^(n2-1)+c2+1;
        id2 = Deg*(loc-1)+[1:Deg];
        i2 = Deg*(j-1)+[1:Deg];

        
        tmp = Mat(id1,id2);
        Mat_tmp = Mat_tmp + sparse(i1'*ones(1,Deg),ones(Deg,1)*i2,tmp,DoFs,DoFs);
        end
    end

end