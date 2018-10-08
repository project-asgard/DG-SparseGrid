function Ix=Index_1D(k,level)

if level==0
    Ix=[1:k];
else
    Ix=[k*2^(level-1)+1:k*2^level];
end

end
