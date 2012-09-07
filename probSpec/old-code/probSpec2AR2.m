function [Lam,Var] = probSpec2AR2(om,lamx,varx)

% function [Lam,Var] = probSpec2AR2(om,lamx,varx)
%
% Converts from the probabilistic spectrogram parameters to AR(2)
% parameters.
%
% INPUTS
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]
%
% OUTPUTS
% Lam = AR(2) dynamical parameter [D,2]
% Var = AR(2) innovations noise, [D,1]

D = length(om);
Lam = zeros(D,2);
Lam(:,1) = 2*lamx.*cos(om);
Lam(:,2) = -lamx.^2;
Var = (1+lamx.^2).*varx;