function [Y,X1,X2,A] = sampleMPAD(Params,Dims)

  % function [Y,X1,X2,A] = sampleMPAD(Params,Dims)
  %
  % Samples T time-steps from the multivariate probabilistic amplitude
  % demodulation. See Turner 2010 Chapter 5 for details of the model.
  %
  % The model either uses carriers which are of the AR(2) type or
  % of the probabilistic filter bank type.
  %
  % y_t = \sum_{d=1}^D c_{d,t} a_{d,t} + \sig_y \epsilon_t
  % \epsilon_t ~ \Norm(0,1)
  %
  % Here c_{d,t} are a set of narrow band carrier which are either
  % second order auto-regressive processes (AR(2)) or probabilistic
  % filter bank processes (pFB): 
  %
  % AR(2)
  % c_{d,t} = \lambda_{1,d} c_{d,t-1} +  ...
  %                    \lambda_{2,d} c_{d,t-2} +\sig_d \eta_{d,t}
  % \eta ~ \Norm(0,1)
  %
  % pFB
  % c_{d,t} = real(z_{d,t})
  % z_{d,t} = \lambda_d \exp(i \omega_d) z_{d,t-1} + ...
  %                     \sig_d (\eta_{d,t}+i\eta_{d,t}
  %
  % The a_{d,t} are a set of positive amplitudes or modulators:
  % a_{d,t} = \log(1+exp(\sum_{k} g_{d,k}  x2_{t,k} + \mu_{d}))
  %
  % The 'transformed envelopes' x2_{t,k} are drawn from Gaussian processes:
  % x2_{1:T,k} ~ \Norm(0,cov_k(t,t'))
  % The covariance function of these Gaussian processes is:
  % cov_k(t,t') = exp(-(t-t')^2/2 len_d^2)
  % where the characteristic time-scale of the Gaussian processes
  % is given by len_d
  %
  % 
  % INPUTS  
  % Params = structure containing
  %   G = Generative weights [D,K] 
  %   Lam1 = Carrier dynamics [D,2] or optionally for pFB carriers [D,1]
  %   Var1 = Carrier noise [1,D]
  %   Len2 = Modulator Length-Scales [K,1]
  %   Var2 = Modulator Noise [1,K]
  %   Mu2 = Modulator means [K,1]
  %   vary = observation noise  
  %   optionally for pFB carriers: om = centre frequencies [D,1]
  %
  % Dims = Structure containing all the dimensionalities
  % D = number of carriers
  % K = number of modulators  
  % T = number of time-steps to sample  
  %
  % OUTPUTS
  % Y = Observations [T,1]
  % X1 = carriers [T,D] (complex valued for pFB)
  % X2 = modulators [T,K]
  % A = amplitudes [T,D]  
    
  

  [D,K,T] = unpackDimsMPAD(Dims);
  
  % Figure out the carrier type
  if min(size(Params.Lam1))==1
    pFB =1; 
  else
    pFB =0;
  end

  % Read out the parameters and sample the carriers
  if pFB==1
    [G,Lam1,Var1,Len2,Var2,Mu2,vary,om] = unpackParamsMPAD(Params);
    
    % Generate the carriers
    [temp,X1] = samplePFB(Lam1,Var1,om,0,T);
  else
    [G,Lam1,Var1,Len2,Var2,Mu2,vary] = unpackParamsMPAD(Params);
    
    % Generate the carriers
    [temp,X1] = sampleAR2FB(Lam1,Var1,0,T);
  end
    

% Generate the transformed amplitudes
X2 = repmat(NaN,[T,K]);

%Blank = '            ';
for k=1:K
%  fprintf(['\rModulator Progress ',num2str(k),'/',num2str(K),Blank])
  X2(:,k) = sampleGPSE(Var2(k),Len2(k),T);
end

%fprintf('\n')

% Amplitudes
A = log(1+exp(X2*G'+ones(T,1)*Mu2'));
  
% Generate the observations
Y = sum(A.*real(X1),2) + randn(T,1)*sqrt(vary);
