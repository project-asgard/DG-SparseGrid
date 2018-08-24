Dim = 4;
Lev = 10;
k = 3;

result = perm_leq(Dim,Lev);
count = 0;
for i = 1:size(result,1)
    tmp = result(i,:);
    count_loc = 1;
    for j = 1:Dim
        count_loc=count_loc*2^max(tmp(j)-1,0);
    end
    count=count+count_loc;
end

FG = (2^Lev*k)^Dim;
SG = k^Dim*count;
[FG SG]