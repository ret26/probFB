function x = sampleMA(phi,T)

% function x = sampleMA(phi,T)
%
% Returns T samples from a moving average process:
%
% x_t = phi_1*e_t + phi_2*e_{t-1}+ ... + phi_tau*e_{t-tau+1}
%
% INPUTS
% phi = moving average parameters size [tau,1]
% T = number of samples
%
% OUTPUTS
% x = samples size [T,1]


tau = length(phi);
noi = randn(T+tau,1);
x = zeros(T,1);
for t=1:T
  x(t) = phi'*noi(t+tau-1:-1:t);
end
