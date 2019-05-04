function [NewInd,NewDof] = RefineElem(Elem,Cel)
    NewElem = Elem+1;
    Cel = [2*Cel,2*Cel+1];
    NewInd = ;
    NewDof = ;
end