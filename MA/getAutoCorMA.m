function autoCor = getAutoCorMA(phi)

% function autoCor = getAutoCorMA(phi)
%
% Return the auto correlation of a moving average process

T = length(phi);

autoCor = zeros(T,1);

for t=1:T
  ind = 1:T-(t-1);
  autoCor(t) = phi(ind)'*phi(ind+(t-1));
end
