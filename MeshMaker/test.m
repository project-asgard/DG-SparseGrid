function test
clc
clear all

%Creat an instance of the class
MM = Mesh_interface

% initial the problem
% Dim, Np, K, Km
MM.init(2, 2, 1, 10)

%find values with keys
keys = [[1 1 0 0 1 0];[1 0 0 0 1 0]];
values = findvs(MM,keys)
values = MM.findvs(keys)

%find a keys with value(only in leaf map)
value = 15
key = findk(MM, value)

key = MM.findk(value)

%refine the problem
%first is elements that need to be added
%second is elements that need to be removed
adds = [15, 12]
dims = [[1,0,0];[0,1,0]];
removes = []
int Nmax = 3;
MM.refine (adds, dims, removes, Nmax)

%print the problem
MM.print
print(MM)

%print the values in the leaf map
MM.leafv
leafv(MM)

%print the values in the all map
MM.allv
allv(MM)

%size of BigMap, number of all keys
MM.size
size(MM)
end
