function [Obj,varargout] = get_Obj_pSTFT_spec(theta,vary,specTar,minVar,limOm,limLam,bet);

% function Obj =
% get_Obj_pSTFT_spec(theta,vary,specTar,minVar,limOm,limLam,bet);
%
% Optional:
%
% function [Obj,dObj] =
% get_Obj_pSTFT_spec(theta,vary,specTar,minVar,limOm,limLam,bet);
%
% For fitting a probabilistic spectrogram model to a signal.  
%
% x_{1,t,d} = lam_d x_{1,t-1,d} +\eta_{1,t,d} \varx_d^{1/2}
% x_{2,t,d} = lam_d x_{2,t-1,d} +\eta_{2,t,d} \varx_d^{1/2}
% \eta_{i,t,d} ~ \Norm(0,1)
% y_t = real(\sum_{d} exp(i om_d)*(x_{1,t,d}+i x_{2,t,d}))

% Fits the parameters (lam_d, varx_d,om_d) such that the spectra are
% matched to the target spectrum specified in specTar.
%
% The cost function is: \sum_t [ log(spec_t) + specTar_t/spec_t ]
% Which can be derived from the fact that the model is a Gaussian
% with a spectrum parameterised using the AR(1)/cosine parameters.
%
% The function parameterises the AR processes using the marginal
% variance of the processes rather than the conditional (marginal
% variance = conditional variance/(1-lam^2)) since the objective is
% more simply behaved in this space and the constraints on the
% parameters are simpler to enforce.
%
% INPUTS
% theta = parameters of the probabilistic spectrogram to be optimised:
%         first D components are the log marginal variances, the next
%         D components are the transformed sinusiod frequencies, the
%         final D parameters the transformed AR lambda
%         parameters. Size is therefore [3*D,1]
% vary = white noise variance (if this is set to zero it can result
%        in numerical problems arising from the division
%        spectAR2./specTar) a sensible setting is: vary=max(specTar)*1e-4
% specTar = target spectrum of the data to be fit, size [N,1]
% minVar = minimum of the marginal variances size [D,1]
% limOm = min/max centre-frequencies e.g. [0,1/2], size [D,2]
% limLam = min/max bandwidths e.g. [0,0.999], size [D,2]
% bet = strength of the gamma shrinkage prior on the marginal
%       variance parameter (set to zero for no pruning)
% 
% OUTPUTS
% Obj = objective
% OPTIONAL OUTPUT:
% dObj = derivative of the objective wrt the parameters, size [3*D,1]
% 
% See Chapter 5 of my thesis (Statistical Models for Natural Sounds
% by R.E.Turner) for more details about AR(2) Filter banks.

D = length(theta)/3;

dVar= exp(theta(1:D));
mVar = minVar+dVar;

om = limOm(:,1)+(limOm(:,2)-limOm(:,1))./(1+exp(-theta(D+1:2*D)));
lam = limLam(:,1)+(limLam(:,2)-limLam(:,1))./(1+exp(-theta(2*D+1:3*D)));

N = length(specTar);
omegas = linspace(0,pi,ceil(N/2));
omegas = [omegas,-omegas([floor(N/2):-1:1])];

% Get the component spectra
spec = ones(1,N)*vary;
for d=1:D

  % for the objective
  const = 1/2*mVar(d)*(1-lam(d)^2);
  alp1 = 1+lam(d).^2-2*lam(d)*cos(omegas-om(d));
  alp2 = 1+lam(d).^2-2*lam(d)*cos(omegas+om(d));  
  spec = spec + const./alp1 + const./alp2; 
    
end


% Likelihood cost function
Obj1 = sum(log(spec))+sum(specTar'./spec);

% Prior cost function
Obj2 = bet*sum(mVar);

% Rescale to get nice units
Obj = (Obj1+Obj2)/N;

%plot(log(specTar),'-k'); hold on; plot(log(spec),'-r')
%keyboard

if nargout>1
  % Derivative of the components wrt parameters

  dObjdtransVar = zeros(D,1);
  dObjdtransOm = zeros(D,1);
  dObjdtransLam = zeros(D,1);
  
  dObjdspec = 1./spec-specTar'./(spec.^2);
  
  for d=1:D
 
    % Derivative wrt transformed marginal variance
    const2 = 1/2*(1-lam(d)^2);
    alp1 = 1+lam(d).^2-2*lam(d)*cos(omegas-om(d));
    alp2 = 1+lam(d).^2-2*lam(d)*cos(omegas+om(d));  
    dspecdtransVar = (mVar(d)-minVar(d))*const2*(1./alp1 + 1./alp2); 
    dObjdtransVar(d) = sum(dObjdspec.*dspecdtransVar); 
    
    % Derivative wrt transformed centre frequency
    const1 = 1/2*mVar(d)*(1-lam(d)^2);
    dspecdom = const1*2*lam(d)*(sin(omegas-om(d))./alp1.^2 ...
				-sin(omegas+om(d))./alp2.^2);
    domdtransOm = (limOm(d,2)-limOm(d,1))*(1/4)./cosh(theta(D+d)/2).^2;
    dObjdtransOm(d) = sum(dObjdspec.*dspecdom)*domdtransOm; 
    
    % Derivative wrt transformed lambda
    const3 = -mVar(d)*lam(d);
    dspecdlam = const3*(1./alp1+1./alp2)+ ...
	     const1*2*((cos(omegas-om(d))-lam(d))./alp1.^2 ...
		       +(cos(omegas+om(d))-lam(d))./alp2.^2);
    dlamdtransLam = (limLam(d,2)-limLam(d,1))*(1/4)./cosh(theta(2*D+d)/2).^2;
    dObjdtransLam(d) = sum(dObjdspec.*dspecdlam)*dlamdtransLam;  
  end
  
  dObj1 = [dObjdtransVar;dObjdtransOm;dObjdtransLam;];

  dObj2 = [bet*dVar;zeros(2*D,1)];
  
  dObj = [dObj1+dObj2]/N;
  varargout{1}= dObj;

end

