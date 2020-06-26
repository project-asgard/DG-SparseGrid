function [Hash,IHash,FineIndex,DPoint] = CoarseMesh(Hash,IHash,Deg,index,FineIndex)


% index = unique(ceil(( MarkedIndex )/Deg));
NewHash = Hash;
NewIHash = IHash;

count_leaf = 1;
NewFineIndex= [];

count = size(IHash,2);
delem = 0;
DPoint = [];
ix_fine = [];

count_del = 1;
for i = 1:size(index,1)
    %     ll = IHash{ceil(index(i)/Deg)}
    ll = IHash{index(i)};
    lev_loc = ll(1);
    cel_loc = ll(2);
    id_loc = ll(3);
    
    key = [lev_loc,cel_loc];
    po = NewHash.(sprintf('i%g_',key));
    NewHash = rmfield(NewHash,sprintf('i%g_',key));
    NewIHash(index(i)-delem) = [];
    DPoint(count_del) = Hash.(sprintf('i%g_',key));
    count_del = count_del +1;
    
%     if lev_loc>1
        ix = find(FineIndex == po);%NewHash.(sprintf('i%g_',key)));
        FineIndex(ix) = [];
        FineIndex(ix:end) = FineIndex(ix:end)-1;
        
        % check about children
        key1 = [lev_loc,2*ceil(cel_loc/2)];
        key2 = [lev_loc,2*ceil(cel_loc/2)+1];
        if isfield(Hash,sprintf('i%g_',key)) == 0 && isfield(Hash,sprintf('i%g_',key)) == 0
            % no children then put mom to fineindex
            key = [lev_loc-1,ceil(cel_loc/2)];
            lr = NewHash.(sprintf('i%g_',key));
            NewFineIndex(count_leaf) = lr;
            count_leaf = count_leaf+1;
        end
%     end
    
    count = index(i)-delem;
    for j = count : size(NewIHash,2)
        ll2 = NewIHash{j};
        key = [ll2(1) ll2(2)];
        NewHash = setfield(NewHash,sprintf('i%g_',key),count);
        count = count+1;
    end
    %     NewHash
    delem = delem + 1;
end

Hash = NewHash;
IHash = NewIHash;


FineIndex=sort([FineIndex,NewFineIndex]);


% FineGrid = EndGrid+1 : Deg*count;