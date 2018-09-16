function xgrid = plotgrid(IHash,Lstart,Lend,yval,option)

% plot grids
count = 0;
for i = 1:size(IHash,2)
    ll = IHash{i};
    h_loc = (Lend-Lstart)*2^(-ll(1)+1);
    count = count+1;
    if ll(1) == 0
        xgrid(count) = Lstart + h_loc/4;
    else
        xgrid(count) = Lstart + h_loc*(ll(2))+h_loc/2;
    end
%     [h_loc Lstart + h_loc*(ll(2))+h_loc/2]
    
end
plot(xgrid,yval*ones(1,count),option)
