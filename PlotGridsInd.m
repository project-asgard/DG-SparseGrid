function PlotGridsInd(Lstart,Lend,IndNZ,VecFlag)

for i = 1:size(IndNZ)
    if i == 1
        plot((Lend+Lstart)/2,0,'ko','MarkerFaceColor','k','MarkerSize',10)
    else
        tmp = IndNZ(i);
        LevT = ceil(log2(tmp));
        CelT = tmp-1-2^(LevT-1);

        hT = (Lend-Lstart)/2^(LevT-1);
        xT(i) = Lstart+CelT*hT+hT/2;
        if VecFlag(tmp) == 3
            plot(xT(i),-LevT,'bo','MarkerFaceColor','b','MarkerSize',10);
        else
            plot(xT(i),-LevT,'ro','MarkerFaceColor','r','MarkerSize',10);
        end
        
    end
    
end