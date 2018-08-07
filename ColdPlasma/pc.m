function y=pc(x)
%------------------------------------------------------------
% Code for preconditioner
% compute M\r
% with M as the preconditioner matrix
%------------------------------------------------------------
global domain Nelem ng nc newlccu newlrru newlcc lnvp01 nlrr nlcc 
global newlb ngB Uc Lc ngB01 nc01 unc1 lnv2  pnc ng1

      
y1=zeros(ng,1);
y2=zeros(nc,1);
yc=zeros(nc,1);
y3=zeros(ng,1);
p0=zeros(pnc,1);

lx=zeros(lnv2,1);
y=zeros(ngB01,1);
llx=zeros(ngB,1);
llx(1:ngB01)=x;
for i=1:Nelem
    lx(newlb)=domain(i).B(newlb,:)*x(1:ng);
    ldx   =domain(i).D(newlrru,newlrru)*lx(newlrru);
    ytemp= domain(i).Ur\(domain(i).Lr\([ldx;zeros(lnvp01,1)]));
    y1= y1 + domain(i).B(newlrru,:)'*domain(i).D(newlrru,newlrru)*ytemp(1:nlrr);
    y2= y2 + domain(i).Bc(newlcc,:)'*(domain(i).Krc)'*ytemp;
    yc= yc + domain(i).Bc(newlcc,:)'*([domain(i).D(newlccu,newlccu)*lx(newlccu);llx(ng+i)]);
end


yc= yc- y2;


y2 = Uc\(Lc\(yc(1:nc01)));

y2(nc)=0;




for i=1:Nelem
    ytempc=domain(i).Bc(newlcc,:)*y2;
    ytemp=domain(i).Ur\(domain(i).Lr\(domain(i).Krc*ytempc));
    ytemp= domain(i).B(newlrru,:)'*domain(i).D(newlrru,newlrru)*ytemp(1:nlrr);
    y3   =y3  + domain(i).B(newlccu,:)'*domain(i).D(newlccu,newlccu)*ytempc(1:nlcc)-ytemp;
end

y(1:ng)= y1+y3;

y(ng1:ngB01)=y2(unc1:nc01);


return;
