
function key=GenerateKey4Dv2(I1,I2,I3,I4,Key1dMesh)


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

row = zeros(1,3*dim);
idx1 = [1:dim:3*dim];
idx2 = [2:dim:3*dim];
idx3 = [3:dim:3*dim];
idx4 = [4:dim:3*dim];

use_loop = (n4 <= 8);
i4 = reshape(1:n4,n4,1);
for i1=1:n1
    row(idx1) = tmp_1(i1,:);
    for i2=1:n2
        row(idx2) = tmp_2(i2,:);
        for i3=1:n3
            row(idx3) = tmp_3(i3,:);

            if (use_loop),
             for i4=1:n4
                row(idx4) = tmp_4(i4,:);
                key(count,:) = row;

                % key(count,[1:dim:3*dim])=tmp_1(i1,:);
                % key(count,[2:dim:3*dim])=tmp_2(i2,:);
                % key(count,[3:dim:3*dim])=tmp_3(i3,:);
                % key(count,[4:dim:3*dim])=tmp_4(i4,:);
                count=count+1;
             end
            else
             mat = repmat(row, [n4,1]);
             mat(i4,idx4) = tmp_4(i4,:);
             key(count:(count+n4-1),1:(3*dim)) = mat(1:n4,1:(3*dim));
             count = count + n4;
            end;
        end
    end
end

key(:,3*dim-dim+1:3*dim)=key(:,3*dim-dim+1:3*dim)-1;
end
