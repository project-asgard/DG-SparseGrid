function FMWT_COMP = OperatorTwoScale(pde,dimIdx,deg,lev)

%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************
% global OperatorTwoScale_method;
idebug = 0;

% imethod_default = 'nonwavelet';
% imethod_default = 'wavelet';
imethod_default = 'wavelet2';

imethod = imethod_default;

% if exist('OperatorTwoScale_method','var'),
%         imethod = OperatorTwoScale_method;
% end;
if (idebug >= 1),
    disp(sprintf('OperatorTwoScale:maxDeg=%d,maxLev=%d,imethod=%s', ...
                                   deg,   lev,   imethod ));
end;

if (strcmp(imethod, 'nonwavelet')),
    FMWT_COMP = OperatorTwoScale_nonwavelet(pde,deg,2^lev);
elseif (strcmp(imethod,'wavelet2')),
    FMWT_COMP = OperatorTwoScale_wavelet2(pde,dimIdx,deg,lev);
else
    FMWT_COMP = OperatorTwoScale_wavelet(pde,dimIdx,deg,2^lev);
end;

end
