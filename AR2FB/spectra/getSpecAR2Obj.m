function [Obj,varargout]=getSpecAR2Obj(theta,vary,specTar);

% function Obj = getSpecAR2Obj(theta,specTar);
% Optional:
% function [Obj,dObj] = getSpecAR2Obj(theta,specTar);
%
% For fitting an AR(2) filter bank to a signal.
% Fits the AR(2) parameters such that the spectra are matched to
% the target spectrum specified in specTar.
%
% The cost function is:
% \sum_d [ log(specAR2_d) + specTar_d/specAR2_d ]
% Which can be derived from the fact that the model is a Gaussian
% with a spectrum parameterised using the AR(2) parameters.
%
% Each AR(2) process is:
% x_{d,t} = \lambda_1 x_{d,t-1} + \lambda_2 x_{d,t-2} + \sigma_d \epsilon
% \epsilon ~ Norm(0,1)
%
% INPUTS
% theta = parameters of the AR2FB to be optimised: first D components
%         are the log-variances (\log \sigma_d, see above), the next D
%         components are the AR(1) dynamical parameters
%         (\lambda_{d,1}) the final D parameters are the
%         \lambda_{d,2}s. Size is therefore [3*D,1]
% vary = white noise variance (if this is set to zero it can result
%        in numerical problems arising from the division
%        spectAR2./specTar) a sensible setting is: vary=max(specTar)*1e-4
% specTar = target spectrum of the data to be fit, size [N,1]
% 
% OUTPUTS
% Obj = objective
% dObj = derivative of the objective wrt the AR(2) parameters, size [3*D,1]
% 
% See Chapter 5 of my thesis (Statistical Models for Natural Sounds
% by R.E.Turner) for more details about AR(2) Filter banks.

D = length(theta)/3;  
Var = exp(theta(1:D));
Lam1 = theta(D+1:2*D);
Lam2 = theta(2*D+1:3*D);

N = length(specTar);
Freq = linspace(0,1/2,ceil(N/2));
Freq = [Freq,-Freq([floor(N/2):-1:1])];

% Get the component spectra
A = (1+Lam1.^2+Lam2.^2)*ones(1,N);
B = 2*Lam1.*(Lam2-1)*cos(2*pi*Freq);
C = -2*Lam2*cos(4*pi*Freq);
specs = Var*ones(1,N)./(A+B+C);  % component spectra
spec = sum(specs,1)+vary; % total spectra

%hold on
%plot(specTar,'-k')
%plot(spec,'-r')
%keyboard
% Objective

Obj = sum(log(spec))+sum(specTar'./spec);
Obj = Obj/N;
%plot(log(specTar),'-k'); hold on; plot(log(spec),'-r')
%keyboard

if nargout>1
  % Derivative of the components wrt parameters

  dspecdtheta = zeros(3*D,N);

  dspecdtheta(1:D,:) = specs;
  
  dspecdtheta(D+1:2*D,:) = -2*specs.^2./(Var*ones(1,N)).*...
      (Lam1*ones(1,N)+(Lam2-1)*cos(2*pi*Freq));
  
  dspecdtheta(2*D+1:3*D,:) = -2*specs.^2./(Var*ones(1,N)).* ...
      (Lam2*ones(1,N)+Lam1*cos(2*pi*Freq)-ones(D,1)*cos(4*pi*Freq));

  dObj = sum(ones([3*D,1])*(1./spec.*(1-specTar'./spec)).* ...
	     dspecdtheta,2);
  
  dObj = dObj/N;
  varargout{1}= dObj;
end