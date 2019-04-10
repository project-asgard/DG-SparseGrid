function y=Fc(x)
% compute A*x
global A_encode

y=zeros(ngB01,1);
y0 =zeros(ngB,1);
yy=zeros(lnv2,1);
yyp=zeros(lnv2,1);
w =zeros(lnv2,1);
xx=zeros(ngB,1);
xx(1:ngB01,1)=x;

for i=1:Nelem
    ipp=i+ng;
    w(newlb)= domain(i).B(newlb,:)*x(1:ng);
    Ui=domain(i).Ui;
    Li=domain(i).Li;
    yy(newlb)= (domain(i).Kbb-domain(i).Kbi*(Ui\(Li\ ...
        domain(i).Kbi')))*w(newlb);
    yyp(newlb)=domain(i).BB0(1,newlb)'*xx(ipp);
    y0(1:ng) = y0(1:ng) +  domain(i).B(newlb,:)'*(yy(newlb)+yyp(newlb));
    y0(ipp)=y0(ipp)+domain(i).BB0(1,newlb)*w(newlb);
end
y=y0(1:ngB01);
return

