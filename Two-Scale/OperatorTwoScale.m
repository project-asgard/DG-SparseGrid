function FMWT_COMP = OperatorTwoScale(maxDeg,maxLev,imethod_in)
%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************

%imethod_default = 'nonwavelet';
imethod_default = 'wavelet';
imethod = imethod_default;
if (nargin >= 3),
        imethod = imethod_in;
end;

if (strcmp(imethod, 'nonwavelet')),
        FMWT_COMP = OperatorTwoScale_nonwavelet(maxDeg,maxLev);
else
        FMWT_COMP = OperatorTwoScale_wavelet(maxDeg,maxLev);
end;

end
