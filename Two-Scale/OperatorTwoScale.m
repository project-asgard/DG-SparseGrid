function FMWT_COMP = OperatorTwoScale(maxDeg,maxLev)
%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************
global OperatorTwoScale_method;
idebug = 1;

imethod_default = 'nonwavelet';
%imethod_default = 'wavelet';

imethod = imethod_default;

if exist('OperatorTwoScale_method','var'),
        imethod = OperatorTwoScale_method;
end;
if (idebug >= 1),
    disp(sprintf('OperatorTwoScale:maxDeg=%d,maxLev=%d,imethod=%s', ...
                                   maxDeg,   maxLev,   imethod ));
end;

if (strcmp(imethod, 'nonwavelet')),
        FMWT_COMP = OperatorTwoScale_nonwavelet(maxDeg,maxLev);
else
        FMWT_COMP = OperatorTwoScale_wavelet(maxDeg,maxLev);
end;

end
