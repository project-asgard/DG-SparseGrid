function sum = Recur_sum(I)
    if I == 0
        sum = 0;
    else
    sum = I+Recur_sum(I-1);
    end
end