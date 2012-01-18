function [Y,X] = sampleAR2FB(Lam,Var,vary,T);
  
  % function [Y,X] = sampleAR2FB(Lam,Var,vary,T);
  %
  % Samples T time-steps from an AR(2) Filter Bank process.  Uses the
  % spectrum of the AR(2) process and filters white Gaussian
  % noise. See sampleAR2FBSlow.m for the auto-regressive version.
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
  % see test_sampleAR2FB.m for the unit tests

  D = length(Var);  
  X = repmat(0,[T,D]);
  Y = repmat(NaN,[T,1]);
    
  % Sample from the AR process by filtering coloured noise
  
  RngFreqs = [0,1/2];
    
  for d=1:D

    % Get spectrum centre frequency
    [Freqs,spec,fMAX,SpecMAX,dF1,dF2] = getSpecAR2(Lam(d,:),Var(d), ...
							    1,RngFreqs);

    tau = 5*ceil(1/fMAX); % if tau=0 the first and last samples
                          % will be correlated - the fft induces
                          % circular correlations
			  
    Tx = 2^ceil(log2(T+tau)); % make a power of two so fft is fast
    
    NumFreqs = Tx/2+1;
    
    % Get spectrum centre frequency
    [Freqs,spec,fMAX,SpecMAX,dF1,dF2] = getSpecAR2(Lam(d,:),Var(d), ...
							    NumFreqs,RngFreqs);
    
    
    spec = [spec,spec(end-1:-1:2)]';
    
    % If you want to handle possibly odd length signals
    % if mod(T,2)==0 
    %   % even
    %   spec = [spec,spec(end-1:-1:2)]';
    % else 
    %   % odd
    %   spec = [spec,spec(end:-1:2)]';
    % end
    
    xCur = ifft(sqrt(spec).*fft(randn(Tx,1)));
    X(:,d) = xCur(1:T);
  end
      
  % Generate the observations from the xs
  Y = sum(X,2) + randn(T,1)*sqrt(vary);