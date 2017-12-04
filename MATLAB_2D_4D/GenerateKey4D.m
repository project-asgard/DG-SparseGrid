function key=GenerateKey4D(I1,I2,I3,I4,Key1dMesh)

tmp_1=Key1dMesh(I1,:);
tmp_2=Key1dMesh(I2,:);
tmp_3=Key1dMesh(I3,:);
tmp_4=Key1dMesh(I4,:);

n1=size(tmp_1,1);
n2=size(tmp_2,1);
n3=size(tmp_3,1);
n4=size(tmp_4,1);
dim=4;
ntol=n1*n2*n3*n4;
key=zeros(ntol,3*dim);
count=1;
for i1=1:n1
    for i2=1:n2
        for i3=1:n3
            for i4=1:n4
                key(count,[1:dim:3*dim])=tmp_1(i1,:);
                key(count,[2:dim:3*dim])=tmp_2(i2,:);
                key(count,[3:dim:3*dim])=tmp_3(i3,:);
                key(count,[4:dim:3*dim])=tmp_4(i4,:);
                count=count+1;
            end
        end
    end
end

key(:,3*dim-dim+1:3*dim)=key(:,3*dim-dim+1:3*dim)-1;
end
