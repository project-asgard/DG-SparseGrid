function MarkedIndex = mark(sol,Index,Deg,epsilon,type)

if type == 'refine'
    MarkedIndex = find(abs(sol(Index))>epsilon);
    
elseif type == 'coarse'
    MarkedIndex = find(abs(sol(Index))<epsilon);
    MarkedIndex = Index(MarkedIndex);
end

MarkedIndex = unique(ceil(MarkedIndex/Deg));

