function [Hash,IHash,FineIndex] = RefineMesh(Hash,IHash,Deg,index,FineIndex)


count = size(IHash,2);
count_Leaf = 1;
NewFineIndex = [];

for i = 1:size(index,1)
    %     ll = IHash{ceil(index(i)/Deg)}
    ll = IHash{index(i)};
    lev_loc = ll(1);
    cel_loc = ll(2);
    if lev_loc < 8%10
        
        if lev_loc>0
            key = [lev_loc+1,2*cel_loc];
        else
            key = [lev_loc+1,0];
        end
        if isfield(Hash,sprintf('i%g_',key)) == 0
            %         'adding 1'
            count = count+1;
            Hash = setfield(Hash,sprintf('i%g_',key),count);
            IHash{count} = [key,2^(lev_loc)+2*cel_loc+1];
            
            ix = find(FineIndex == index(i));
            
            if size(ix,1)>0
                FineIndex(ix) = [];
            end
            NewFineIndex(count_Leaf) = count;
            count_Leaf = count_Leaf+1;
            %         end
        end
        
        
        if lev_loc>0
            key = [lev_loc+1,2*cel_loc+1];
        else
            key = [lev_loc+1,0];
        end
        if isfield(Hash,sprintf('i%g_',key)) == 0
            %         'adding 2'
            count = count+1;
            Hash = setfield(Hash,sprintf('i%g_',key),count);
            IHash{count} = [key,2^(lev_loc)+2*cel_loc+2];
            
            ix = find(FineIndex == index(i));
            if size(ix,1)>0
                FineIndex(ix) = [];
                
            end

            NewFineIndex(count_Leaf) = count;
            count_Leaf = count_Leaf+1;

        end
    end
    
end

FineIndex = sort(unique([FineIndex,NewFineIndex]));
% FineGrid = EndGrid+1 : Deg*count;