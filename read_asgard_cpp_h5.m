function data = read_asgard_cpp_h5(filename)

hinfo = hdf5info(filename);
dset = hdf5read(hinfo.GroupHierarchy.Datasets(1));

data = dset(:,end);

n = numel(data);
n2 = sqrt(n);
data2d = reshape(data,n2,n2);
figure
contourf(data2d)
figure
semilogy(data2d(end,:))
hold on
semilogy(data2d(1,:))

end
