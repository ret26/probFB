function [fftCov,varargout] = getGPSEisoSpec(sig2,len,T)
%getGPSMSpec gets spectrum of a spectral mixture kernel
%  Inputs:
%    var: variance
%    len: lengthscale
%    T: length of data
%  Outputs:
%    fftCov: spectrum of the kernel, evaluated at [0:T/2 -T/2:1]
%    dfftCov: optional, derivatives of the kernel wrt hypers
%
%  k(f) = sqrt(2*pi)*sig2*len*(exp(-2*pi^2*len^2*f.^2);
%
% See also getGPSESpec.m, getGPSMSpec.m
% 
% Written by Richard Turner and Thang Bui
% Last modified: 10/2014

f = [0:ceil(T/2) -floor(T/2)+1:1:-1]'/T;

fftCov = sqrt(2*pi)*sig2*len*exp(-2*pi^2*len^2*f.^2);
% derivatives with respect to log of vars, lens and freqs
if nargout>1
    dfftCov = zeros(T,2);
    dfftCov(:,1) = fftCov; % d/dlogvarx_k
    M = -4*pi^2*len^2 * fftCov;
    dfftCov(:,2) = fftCov + f.^2.*M;
    varargout{1} = dfftCov;
end

end