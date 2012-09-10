function [fmax,df,varMa] = probSpec2freq(om,lamx,varx)

% function [fmax,df,varMa] = probSpec2freq(om,lamx,varx)
%
% Finds the spectral parameters of a probabilistic spectrogram.
%
% INPUTS
% lamx = dynamical AR parameters [D,1]
% varx = dynamical noise parameters [D,1]
% om = mean frequencies of the sinusoids [D,1]
%
% OUTPUTS
% fmax = centre frequencies, size [D,1]
% df = bandwidths, size [D,1]
% varMa = marginal variances, size [D,1]

fmax = om/(2*pi);
df = acos(2-(lamx.^2+1)./(2*lamx))/(2*pi);
varMa = varx./(1-lamx.^2);


%[Lam,Var] = probSpec2AR2(om,lamx,varx);
%[fmax,df, varMa] = AR22freq(Lam,Var);