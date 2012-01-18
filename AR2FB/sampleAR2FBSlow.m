function [Y,X] = sampleAR2FBSlow(Lam,Var,vary,T);
  
  % function [Y,X] = sampleAR2FBSlow(Lam,Var,vary,T);
  %
  % Samples T time-steps from an AR(2) Filter Bank process.  Uses the
  % Auto-regressive form of the model - this involves looping over
  % time steps and is therefore a slower method of sampling than using
  % the spectrum of the AR(2) process and filtering white Gaussian
  % noise. See sampleAR2FB.m for the fast version.
  %
  % Inputs
  % Lam = initial dynamical parameters [D,2]
  % Var = initial dynamic noise variance [D,1]
  % vary = initial observation noise 
  % T = number of time-steps to sample
  %
  % Outputs
  % Y = observations [T,1] 
  % X = latent variables [T,D]
  %
  % see Turner, 2010 Chapter 5 for details of the AR(2) filter bank 
  % see test_sampleAR2FBSlow.m for the unit tests

  D = length(Var);  
  X = repmat(0,[T,D]);
  Y = repmat(NaN,[T,1]);
    
  % Sample the first two xs from the stationary distribution of the
  % AR(2) process
  for d=1:D
    autoCor = getAutoCorARTau(Lam(d,:)',Var(d),2);
    
    covX = [autoCor';autoCor(2:-1:1)']; 

    X(1:2,d) = chol(covX)'*randn(2,1);
  end

  % Sample the rest of the xs using the auto-regressive formulation
  for t = 3:T
    X(t,:) = Lam(:,1)'.*X(t-1,:) + Lam(:,2)'.*X(t-2,:) ...
                                 + sqrt(Var)'.*randn(1,D);
  end
      
  % Generate the observations from the xs
  Y = sum(X,2) + randn(T,1)*sqrt(vary);
  
