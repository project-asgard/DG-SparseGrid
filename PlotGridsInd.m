function PlotGridsInd(Lstart,Lend,IndNZ,VecFlag)

% for i = 1:size(IndNZ)
%     if i == 1
%         plot((Lend+Lstart)/2,0,'ko','MarkerFaceColor','k','MarkerSize',10)
%     else
%         tmp = IndNZ(i);
%         LevT = ceil(log2(tmp));
%         CelT = tmp-1-2^(LevT-1);
% 
%         hT = (Lend-Lstart)/2^(LevT-1);
%         xT(i) = Lstart+CelT*hT+hT/2;
%         if VecFlag(tmp) == 3
%             plot(xT(i),-LevT,'bo','MarkerFaceColor','b','MarkerSize',10);
%         elseif VecFlag(tmp) == 2
%             plot(xT(i),-LevT,'go','MarkerFaceColor','g','MarkerSize',10);
%         elseif VecFlag(tmp) == 1
%             plot(xT(i),-LevT,'ro','MarkerFaceColor','r','MarkerSize',10);
%         end
%         
%     end
%     
% end

for i = 1:size(IndNZ)
    if i == 1
        plot((Lend+Lstart)/2,-.1,'k.','MarkerFaceColor','k')
    else
        tmp = IndNZ(i);
        LevT = ceil(log2(tmp));
        CelT = tmp-1-2^(LevT-1);

        hT = (Lend-Lstart)/2^(LevT-1);
        xT(i) = Lstart+CelT*hT+hT/2;
        if VecFlag(tmp) == 3
            plot(xT(i),-.1,'b.','MarkerFaceColor','b');
        else
            plot(xT(i),-.1,'r.','MarkerFaceColor','r');
        end
        
    end
    
end