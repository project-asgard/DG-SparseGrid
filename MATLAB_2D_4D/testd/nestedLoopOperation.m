function counters = nestedLoopOperation(counters, lengths, level)
    if level == length(counters)+1
        %do something
        disp(counters)
        
    else
        for i = 0:lengths(level)
          counters(level) = i;
          counters = nestedLoopOperation(counters, lengths, level + 1);
        end
    end
end