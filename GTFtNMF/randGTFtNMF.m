function [y,V,H,Z] = randGTFtNMF(W,lamv,varv,omv,vary,lenx,mux,varx,tol,T);

% function [y,V,H,Z] = randGTFtNMF(W,lamv,varv,omv,vary,lenx,mux,varx,tol,T);
% 
% Draws variables from the GTF-tNMF model using the hierarchical
% filter bank formulation 
%
% y_t = sum_d Re(v_{d,t}) + sigy*\eta_t
% v_{d,t} = amp_{d,t}/amp_{d,t-1}*lamv_d*exp(i omv_d)*v_{d,t-1} 
%                              + amp_{d,t} \sigma_d \epsilon_{d,t}
% amp_{d,t} = sqrt(H*W/2)
% W_{k,d} = exp(logW_{k,d})/sum_d exp(logW_{k,d})
%
% Each of the temporal basis functions in NMF is given by: 
% H_{t,k} = exp(X_{t,k})  where X{1:T,k} ~ mux_k + GP(lenx_k,varx_k)
%
% The basis function representation of the GP used for the temporal
% basis functions is employed rather than the function-space view
% so that:
%
% bas_k \propto sqrt(varx_k)*exp(-1/(lenx_k^2) t^2)
% X_{1:T,k} = convolve(Z_{1:K,k},bas_k) 
% 
%
% INPUTS
% W = spectral weights [K,D]
% lamv = dynamical AR parameters [D,1]
% varv = dynamical noise parameters [D,1]
% omv = mean frequencies of the sinusoids [D,1]
% vary = observation noise [T,1]
% lenx = squared exponential length-scale [1,K]
% mux = steady state mean of the log temporal priors [1,K]
% varx = steady state variance of the log temporal priors [1,K]
% tol = tolerance parameter which sets the extent of the Gaussian
%       basis functions
%
% OUTPUTS
% y = signal [T,1]
% V = filter bank coefficients [T,D]
% H = temporal basis functions from the NMF component [T,D]
% Z = matrix of temporal basis functions [T,K]


% read out the dimensions of the model  
D = length(lamv);
K = length(lenx);
  
% get the temporal basis functions
[temp,H] = randtnmf(W,lenx,mux,varx,T);

% produce the envelopes
A = (1/2*H*W).^0.5;

% get the filter banks coefficients
[temp,Z] = samplePFB(lamv,varv,omv,0,T);
V = Z.*A;

% produce the signal
y = sum(real(V),2)+sqrt(vary).*randn(T,1);

