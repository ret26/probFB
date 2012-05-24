function X = AR2FBFFT(y,Lam,Var,vary,varargin)
  
  % function X = AR2FBFFT(y,Lam,Var,vary)
  %
  % Second order auto-regressive filter bank (see Turner 2010,
  % Chapter 5 for details).
  %
  % Implemented using the FFT - fast but only returns the means and
  % cannot handle non-stationary noise.
  %
  % INPUTS
  % y = data [T,1]
  % Lam = dynamical parameters [D,2]
  % Var = dynamic noise variance [D,1]
  % vary = observation noise 
  %  
  % OUTPUTS
  % X = AR2 process mean values, size [D,T]
  %

T = length(y);
  
RngFreqs = [0,1/2];
D = length(Var);


% DETERMINE PADDING BASED ON LONGEST TIME-CONSTANTS OF Xs
tau = 0;
tol = 3;
for d=1:D

  % Get spectrum centre frequency
  [Freqs,spec,fMAX,SpecMAX,dF1,dF2] = getSpecAR2(Lam(d,:),Var(d), ...
							    1,RngFreqs);
  tau = max([0,tol*ceil(1/fMAX)]); % if tau=0 the first and last samples
				 % will be correlated - the fft induces
				 % circular correlations
end

Tx = 2^ceil(log2(T+tau)); % make a power of two so fft is fast
NumFreqs = Tx/2+1;

% GET THE SPECTRA OF THE COMPONENTS
specX = zeros(Tx,D);
for d=1:D
    
  % Get spectrum centre frequency
  [Freqs,spec,fMAX,SpecMAX,dF1,dF2] = getSpecAR2(Lam(d,:),Var(d), ...
							  NumFreqs,RngFreqs);
    
  specX(:,d) = [spec,spec(end-1:-1:2)]';
 
end

% USE THE SPECTRA TO FILTER THE SIGNAL

% edges have to be handled heuristically as they are effectively
% non-stationary

dT = Tx-T;

% linearly interpolate the extra chunk
yExtra = linspace(y(end),y(1),dT)';

% sigmoid interpolation
%len = dT/10;
%yExtra = y(end)+(y(1)-y(end))./(1+exp(-([1:dT]'-dT/2)/len));
%plot(yExtra)
%keyboard

yFFT = fft([y(:);yExtra]);
specY = sum(specX,2)+vary;
X = zeros(T,D);

for d=1:D

    curSpec = specX(:,d)./specY;
  
    xCur = ifft(curSpec.*yFFT);
    X(:,d) = xCur(1:T);
end

X = X'; % for agreement with old method