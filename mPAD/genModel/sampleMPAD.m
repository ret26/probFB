function [Y,X1,X2,A] = sampleMPAD(Params,Dims)

  % function [Y,X1,X2,A] = sampleMPAD(Params,Dims)
  %
  % Samples T time-steps from the multivariate
  % probabilistic amplitude demodulation. See Turner 2010 Chapter 5
  % for details of the model.
  %
  % INPUTS  
  % Params = structure containing
  %   G = Generative weights [D,K] 
  %   Lam1 = Carrier dynamics [D,2]
  %   Var1 = Carrier noise [1,D]
  %   Len2 = Modulator Length-Scales [K,1]
  %   Var2 = Modulator Noise [1,K]
  %   Mu2 = Modulator means [K,1]
  %   vary = observation noise  
  %
  % Dims = Structure containing all the dimensionalities
  % D = number of carriers
  % K = number of modulators  
  % T = number of time-steps to sample  
  %
  % OUTPUTS
  % Y = Observations [T,1]
  % X1 = carriers [T,D]
  % X2 = modulators [T,K]
  % A = amplitudes [T,D]  
    
[G,Lam1,Var1,Len2,Var2,Mu2,vary] = unpackParamsMPAD(Params);
[D,K,T] = unpackDimsMPAD(Dims);

% Generate the carriers
[temp,X1] = sampleAR2FB(Lam1,Var1,0,T);

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
Y = sum(A.*X1,2) + randn(T,1)*sqrt(vary);
