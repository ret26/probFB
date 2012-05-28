function S = probSTFT_FFT(y,lamx,varx,om,vary)
  
  % function S = probSTFT_FFT(y,lamx,varx,om,vary)
  %
  % Probabilistic STFT (see Turner 2010, Chapter 5 for details)
  %
  % Implemented using the FFT - fast but only returns the means and
  % cannot handle non-stationary noise.
  %
  % INPUTS
  % lamx = dynamical AR parameters [D,1]
  % varx = dynamical noise parameters [D,1]
  % om = mean frequencies of the sinusoids [D,1]
  % vary = oberservation noise
  % y = Data, size [T,1] 
  %
  %  
  % OUTPUTS
  % S = probabilistic STFT process mean values (these are complex), size [D,T]
  %
  
S = probFB_FFT(y,lamx,varx,om,vary);

[D,T] = size(S);


for d=1:D
  S(d,:) = exp(-i*om(d)*[1:T]).*S(d,:);
end

