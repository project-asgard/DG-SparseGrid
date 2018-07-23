function num = Start_LevCell(Ix,Iy)
    num = 0;
    for sum = 0:Ix+Iy-1
        for Jx = sum:-1:0
            tmp = Num4Cell(Jx,sum-Jx);
            num = num+tmp;
%             [sum Jx tmp num ]
        end
    end
    for Jx = Ix+Iy:-1:Ix+1
        tmp = Num4Cell(Jx,Ix+Iy-Jx);
        num = num+tmp;
    end
%     if Ix == 0 && Iy == 0
%         num = 0;
%     else
%         num = Num4Cell(Ix-1,Iy-1);
%     end
end

