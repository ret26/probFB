function [om,lamx,varx] = AR22probSpec(Lam,Var)

% function [om,lamx,varx] = AR22probSpec(Lam,Var) 
%
% Converts from the AR(2) parameters to the probabilistic spectrogram
% parameters.
%
% INPUTS
% Lam = AR(2) dynamical parameter [D,2]
% Var = AR(2) innovations noise, [D,1]
%
% OUTPUTS
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]

lamx = sqrt(-Lam(:,2));
om = acos(Lam(:,1)./(2*lamx));
varx = Var./(1+lamx.^2);