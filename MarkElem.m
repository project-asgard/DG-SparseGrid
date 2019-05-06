function ind = MarkElem(f,Type,ETA)
%     tmp = reshape(f,Deg,size(f,1)/Deg);
    
    if strcmp(Type,'refine') == 1
        ind = find(sqrt(sum(f.^2))>ETA);
        
    elseif strcmp(Type,'coarse') == 1
        ind = find(sqrt(sum(f.^2))<ETA);
        
    end
        

end