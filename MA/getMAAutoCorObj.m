function [Obj,dObj] = getMAAutoCorObj(phi,autoCorTar,smth)

% function [Obj,dObj] = getMAAutoCorObj(phi,autoCorTar,smth)
%
% Computes the squared error between a moving average
% autocorrelation and a target auto-correlation and the derivative
% of that squared error. The Moving Average process is defined as:
%
% x_t = \sum_{t'=1}^T phi_t' * \epsilon_{t-t'}
%

T = length(autoCorTar);

autoCor = zeros(T,1);

for t=1:T
  ind = 1:T-(t-1);
  autoCor(t) = phi(ind)'*phi(ind+(t-1));
end
  
delta = autoCor-autoCorTar;
Obj1 = sum(delta.^2)/2;
dObj1 = zeros(T,1);

for t=1:T
  ind1 = 1:T-(t-1);
  ind2 = 1:t;
  dObj1(t) = delta(ind1)'*phi(ind1+t-1)+delta(ind2)'*phi(t-ind2+1);
%keyboard
end

% smoothness term
dPhi = phi(1:T-1)-phi(2:T);
Obj2 = 1/2*sum(dPhi.^2);

dObj2 = [dPhi(1:T-1);0]-[0;dPhi(1:T-1)];

Obj = Obj1+smth*Obj2;
dObj = dObj1+smth*dObj2;