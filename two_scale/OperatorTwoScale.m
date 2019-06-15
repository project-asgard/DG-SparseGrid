function FMWT_COMP = OperatorTwoScale(deg,lev)

%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************
% global OperatorTwoScale_method;
idebug = 0;

%imethod_default = 'nonwavelet';
%imethod_default = 'wavelet';
imethod_default = 'wavelet2';

imethod = imethod_default;

% NOTE: this was overwriting the choice above so commented out.
% if exist('OperatorTwoScale_method','var'), 
%         imethod = OperatorTwoScale_method;
% end;
if (idebug >= 1),
    disp(sprintf('OperatorTwoScale:maxDeg=%d,maxLev=%d,imethod=%s', ...
                                   deg,   lev,   imethod ));
end;

if (strcmp(imethod, 'nonwavelet')),
    FMWT_COMP = OperatorTwoScale_nonwavelet(deg,lev);
elseif (strcmp(imethod,'wavelet2')),
    FMWT_COMP = OperatorTwoScale_wavelet2(deg,lev);
else
    FMWT_COMP = OperatorTwoScale_wavelet(deg,lev);
end;

end
