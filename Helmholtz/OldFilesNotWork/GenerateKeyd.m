function key=GenerateKeyd(Id,Key1dMesh)

dim=size(Id,1);
n=zeros(dim,1);

ntol=1;
for i=1:dim
    tmp{:,i}=Key1dMesh(Id{i},:);
    n(i)=size(tmp{:,i},1);
    ntol=ntol*n(i);
end
ntol
% key=zeros(ntol,3*dim);
count=1;

for dd=dim:-1:1
    for i=1:n(dd)
        key(count,[dd:dim:3*dim])=tmp{:,dd}(i,:);
        count=count+1;
    end
end

% key(:,3*dim-dim+1:3*dim)=key(:,3*dim-dim+1:3*dim)-1;

end
