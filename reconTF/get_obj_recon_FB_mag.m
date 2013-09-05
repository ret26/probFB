function [obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec,flag,varargin)

% function [obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec)
%
% computes error between a signal's subband envelopes and target
% envelopes, and the derivative of this error. Used for
% reconstruction of phaseless filterbanks.
%
% [obj,dobj] = get_obj_recon_FB_mag(y,aTar,spec,flag,lamy,vary)
% optional regularisation which models signals as AR processes
%
% INPUT
% y =  signal [T,1]
% aTar = target envelopes [T,D]
% spec = subband filter ffts [T,D]
% flag = determines the error measure: 
%        'a' squared error in envelopes
%        'loga' squared error in log-envelopes
% optional inputs:
% lamy = AR dynamical parameter
% vary = AR innovations noise parameter
%
% OUPUT
% obj = objective (squared error in evelopes) 
% dobj = derivative in squared error
%
% see test_get_obj_recon_FB_mag for unit tests

[T,D] = size(aTar);
yFFT = fft(y);

obj = 0;
dobj = zeros(T,1);


if strcmp(flag,'a')
  
  for d=1:D
    zCur = ifft(spec(:,d).*yFFT);
    aCur = abs(zCur);
    
    delta_a = aCur-aTar(:,d);
    obj = obj + 0.5*sum(delta_a.^2);
  
    beta = delta_a.*conj(zCur)./aCur;
    dobj = dobj + real(fft(ifft(beta).*spec(:,d)));
  end
  
elseif strcmp(flag,'loga')

  % seems to perform worse by perceptual quality measures than the above
  offset = 1e-3;
  
  for d=1:D
    zCur = ifft(spec(:,d).*yFFT);
    aCur = abs(zCur);
    
    delta_a = log(aCur+offset)-log(aTar(:,d)+offset);
    obj = obj + 0.5*sum(delta_a.^2);
  
    beta = delta_a.*conj(zCur)./(aCur.*(aCur+offset));
    dobj = dobj + real(fft(ifft(beta).*spec(:,d)));
  end

else
  disp('unknown envelope objective');
  return;
end


if nargin>4
  % optional AR prior on signal
  
  lam = varargin{1};
  vary = varargin{2};
  
  obj2 = 1/(2*vary)*lam^2*y(1)^2+1/(2*vary)*y(T)^2 + ...
	1/(2*vary)*(1+lam^2)*sum(y(2:T-1).^2)-...
	lam/vary*sum(y(1:T-1).*y(2:T));

  dobj2 = zeros(T,1);
  dobj2(1) = lam^2/vary*y(1)-lam/vary*y(2);
  dobj2(2:T-1) = (1+lam^2)/vary*y(2:T-1)-lam/vary*(y(3:T)+y(1:T-2));
  dobj2(T) = y(T)/vary-lam/vary*y(T-1);

  obj = obj + obj2;
  dobj = dobj + dobj2;
end
  
obj = obj/(D*T);
dobj = dobj/(D*T);
