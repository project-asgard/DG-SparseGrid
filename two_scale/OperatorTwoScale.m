function FMWT_COMP = OperatorTwoScale(deg,lev,user_method)

%----------------------------------
% Set-up Two-scale operator       %
%----------------------------------
% Input: Degree: maxDeg
%        Level: Np
% Output: Convert Matrix: FMWT_COMP
%**********************************
% global OperatorTwoScale_method;
idebug = 0;

%%
% user_method options
%
% user_method_default = 'wavelet';    % Original method which requires zeroing out
% user_method_default = 'wavelet2';   % Ed's update which only calculates the non-zeros
% user_method_default = 'nonwavelet'; % The non-wavelet approach 

imethod = 'wavelet2';

if exist('user_method','var') && ~isempty(user_method)
    imethod = user_method;
end

if (idebug >= 1),
    disp(sprintf('OperatorTwoScale:maxDeg=%d,maxLev=%d,imethod=%s', ...
                                   deg,   lev,   imethod ));
end

FMWT_COMP = 0;
if (strcmp(imethod, 'nonwavelet'))
    FMWT_COMP = OperatorTwoScale_nonwavelet(deg,lev);
elseif (strcmp(imethod,'wavelet2'))
    FMWT_COMP = OperatorTwoScale_wavelet2(deg,lev);
elseif (strcmp(imethod,'wavelet'))
    FMWT_COMP = OperatorTwoScale_wavelet(deg,lev);
else
    error("Invalid method choice in OperatorTwoScale");
end

end
