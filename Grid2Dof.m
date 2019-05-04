function id = Grid2Dof(IndGrid,Deg)
    tmp = (repmat(IndGrid,1,Deg)-1)*Deg+[1:Deg];
    tmp = tmp';
    id = tmp(:);
    
end