function phi = MAAutoCor2Phi(autoCor,numIts,smth)

% function phi = MAAutoCor2Phi(autoCor,numIts,smth)
%
% Returns the parameters of a moving average process, phi, that
% best matches an autocorrelation, autoCor, in a squared error
% sense. The order of the MA process is set by the length of the
% auto-correlation vector. 
%
% x_t = phi_1*e_t + phi_2*e_{t-1}+ ... + phi_T*e_{t-T+1}
% T = length(autoCor)
%
% The function gives you the option of introducing some smoothing
% of the autoCorrelation - as measured by the squared error between
% adjacent elements in phi - and the variable smth controls the
% weighting and therefore the degree of smoothness.
%
% Note that the squared error metric, although simple to handle,
% would probably be more correct if replaced by a KL objective
%
% INPUTS
% autoCor = auto-correlation, size [T,1]
% numIts = number of iterations, scalar
% smth = degree of smoothing, scalar
% 
% OUTPUTS
% phi = moving average parameters, size [T,1]

%phi = randn(T,1);
phi = sqrt(abs(autoCor)).*sign(autoCor);
[phi,Obj,in] = minimize(phi,'getMAAutoCorObj',numIts,autoCor,smth);
autoCorMA = getAutoCorMA(phi);

% ambiguity up to a flip
if abs(phi(end))>abs(phi(1))
  phi = phi(end:-1:1);
end
