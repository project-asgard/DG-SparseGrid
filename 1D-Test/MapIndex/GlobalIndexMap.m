function index = GlobalIndexMap(Ix,Iy,Deg)

StartIndex = Start_LevCell(Ix,Iy)+1; 
% if we consider the cell information

EndIndex = StartIndex+Num4Cell(Ix,Iy)-1;
% [StartIndex EndIndex]
index = Deg*(StartIndex-1)+1:Deg*(EndIndex);

end

function sum = Recur_sum(I)
    if I == 0
        sum = 0;
    else
    sum = I+Recur_sum(I-1);
    end
end

