
function res = DFS(k,n,start,res,sol)
    if k<0 || n<0
       return;
    else if k==0 && n==0
            res(end+1,:) = sol;
            return;
        else
            for i=start:n
                sol(end+1,:) = i;
                res = DFS(k-1,n-i,start,res,sol);
                sol(end,:) = [];
            end
        end
    end
end
