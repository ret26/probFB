function X = sampleARTau(lams,varx,T,N)
  
  % function x = sampleARTau(lams,varx,T,N)
  %
  % Draw N realisations from an AR(Tau) process T samples long. The
  % first tau samples are drawn from the stationary distribution of
  % the process. 
  %
  % INPUTS
  % lams = dynamical parameters [1,tau]
  % varx = dynamical noise  
  % T = length of time series
  % N = Number of time series
  %  
  % OUTPUTS
  % X = realisations [T,N]
    
  tau = length(lams);
  
  X = repmat(0,[T,N]);
  
  Lams = repmat(lams',[1,N]);
    
  autoCor = getAutoCorARTau(lams',varx,tau);

  covX = zeros(tau);
  for t=1:tau
    covX = covX + diag(ones(tau-t+1,1),t-1)*autoCor(t); 
  end
  
  X(1:tau,:) = chol(covX)'*randn(tau,N);
  
  for t = tau+1:T
    X(t,:) = sum(Lams.*X(t-1:-1:t-tau,:),1) + sqrt(varx).*randn(1,N); 
  end
    
