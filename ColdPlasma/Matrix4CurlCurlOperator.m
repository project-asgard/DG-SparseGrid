function Mat = Matrix4CurlCurlOperator(Lev,Deg,Lmax)
% Compute matrices for curl curl operator
% Operators:: u_h and v denote trial and test functions
% - Volume Integral
%   DuDv = (u_h',v')
%   DuIv = (u_h',v)
%   IuIv = (u_h,v)
% - Interface Integral
%   AuJv = <{u_h},[v]>_F
%   JuAv = <[u_h],{v}>_F
%   ADuJv = <{u_h'},[v]>_F
%   JuADv = <[u_h],{v'}>_F
%   JuJv = <[u_h],[v]>_F
%   JDuJDv = <[u_h'],[v']>_F
% Note:: AuJv = JuAv'; ADuJv = JuADv';
%--------by Lin Mu, 08-09-2018------------

% compute the trace values
p_1 = legendre(-1,Deg);
p_2 = legendre( 1,Deg);

quad_num = 10;

[quad_x,quad_w] = lgwt(quad_num,-1,1);
p_val = legendre(quad_x,Deg);
Dp_val = dlegendre(quad_x,Deg);
Dp_1 = dlegendre(-1,Deg);
Dp_2 = dlegendre(1,Deg);


%---------------------------
% Define dofs
%---------------------------
nx = 2^(Lev);
hx = Lmax/nx;
dofs = Deg*nx;

% - Volume Integral
% Matrix for (u_h',v')
val = Dp_val'*(quad_w*ones(1,Deg).*Dp_val)*1/hx/hx*2;
Ac = repmat({val},nx,1);
DuDv = blkdiag(Ac{:});

% Matrix for (u_h',v)
val = p_val'*(quad_w*ones(1,Deg).*Dp_val)/hx;
Ac = repmat({val},nx,1);
DuIv = blkdiag(Ac{:});

% Matrix for (u_h,v) = eye
IuIv = speye(dofs);

% - Interface Integral
% Matrix for <{u_h},[v]>
Amd  = -p_1'*p_1/2+p_2'*p_2/2;
Asub = -p_1'*p_2/2;
Asup =  p_2'*p_1/2;
AuJv = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);
% Correction for boundary
AuJv(1:Deg,1:Deg) = AuJv(1:Deg,1:Deg)+(-p_1'*p_1)/hx/2;
AuJv(end-Deg+[1:Deg],end-Deg+[1:Deg]) = AuJv(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
    (p_2'*p_2)/hx/2;

% Matrix for <[u_h],{v}>
Amd  = p_1'*(-p_1)/2+p_2'*p_2/2;
Asub =  p_1'*p_2/2;
Asup =  p_2'*(-p_1)/2;
JuAv = 1/hx*blktridiag([Amd],[Asub],[Asup],nx);
% Correction for boundary
JuAv(1:Deg,1:Deg) = JuAv(1:Deg,1:Deg)+(p_1'*(-p_1))/hx/2;
JuAv(end-Deg+[1:Deg],end-Deg+[1:Deg]) = JuAv(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
                (p_2'*p_2)/hx/2;
            
% Matrix for <{u_h'},[v]>ds
Amd = -p_1'*Dp_1/2+p_2'*Dp_2/2;
Asub =-p_1'*Dp_2/2;
Asup = p_2'*Dp_1/2;
ADuJv = 1/hx*blktridiag([Amd],[Asub],[Asup],nx)*2^(Lev+1);
% Correction for boundary
ADuJv(1:Deg,1:Deg) = ADuJv(1:Deg,1:Deg)+(-p_1'*Dp_1)/hx*2^(Lev+1)/2;
ADuJv(end-Deg+[1:Deg],end-Deg+[1:Deg]) = ADuJv(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
    (p_2'*Dp_2)/hx*2^(Lev+1)/2;

% Matrix for <[u_h],{v'}>ds
Amd = Dp_1'*(-p_1)/2+Dp_2'*p_2/2;
Asub = Dp_1'*p_2/2;
Asup = Dp_2'*(-p_1)/2;
JuADv = 1/hx*blktridiag([Amd],[Asub],[Asup],nx)*2^(Lev+1);
% Correction for boundary
JuADv(1:Deg,1:Deg) = JuADv(1:Deg,1:Deg)+(Dp_1'*(-p_1))/hx*2^(Lev+1)/2;
JuADv(end-Deg+[1:Deg],end-Deg+[1:Deg]) = JuADv(end-Deg+[1:Deg],end-Deg+[1:Deg])+...
                (Dp_2'*p_2)/hx*2^(Lev+1)/2;

% Matrix for <[u_h],[v]>/hx
Amd  =  p_1'*(p_1)+p_2'*p_2;
Asub = -p_1'*p_2;
Asup =  p_2'*(-p_1);
JuJv = 1/hx*blktridiag([Amd],[Asub],[Asup],nx)*2/hx;

% Matrix for <[u_h'],[v']>
Amd  =  Dp_1'*(Dp_1)+Dp_2'*Dp_2;
Asub = -Dp_1'*Dp_2;
Asup =  Dp_2'*(-Dp_1);
JDuJDv = 1/hx*blktridiag([Amd],[Asub],[Asup],nx)*2/hx;

% Matrix for <[u_h],[v']>
% JuJDv = 

% convert to multiwavelet basis
FMWT_COMP = OperatorTwoScale(Deg,2^Lev);

DuDv = FMWT_COMP*DuDv*FMWT_COMP';
DuIv = FMWT_COMP*DuIv*FMWT_COMP';
AuJv = FMWT_COMP*AuJv*FMWT_COMP';
JuAv = FMWT_COMP*JuAv*FMWT_COMP';
ADuJv = FMWT_COMP*ADuJv*FMWT_COMP';
JuADv = FMWT_COMP*JuADv*FMWT_COMP';
JuJv = FMWT_COMP*JuJv*FMWT_COMP';
JDuJDv = FMWT_COMP*JDuJDv*FMWT_COMP';

% DuDv(find(abs(DuDv)<1e-5))=0;
% DuIv(find(abs(DuIv)<1e-5))=0;
% AuJv(find(abs(AuJv)<1e-5))=0;
% JuAv(find(abs(JuAv)<1e-5))=0;
% ADuJv(find(abs(ADuJv)<1e-5))=0;
% JuADv(find(abs(JuADv)<1e-5))=0;
% JDuJDv(find(abs(JDuJDv)<1e-5))=0;


Mat = struct('DuDv',DuDv,'DuIv',DuIv,'IuIv',IuIv,'AuJv',AuJv,'JuAv',JuAv,...
    'ADuJv',ADuJv,'JuADv',JuADv,'JuJv',JuJv,'JDuJDv',JDuJDv);