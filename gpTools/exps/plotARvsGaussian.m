% load signal
File = '74 - Sentences'; % Name of file to load
fs = 16000; % sampling rate of file
%RngLimTrain = round([fs*1/2+1,2.1*fs]);  % one sentence
RngLimTrain = round([fs*1/2+1,4.2*fs]);  % two sentences
%RngLimTrain = round([fs*1/2+1,6.5*fs]);  % three sentences

%RngLimTest = round([fs*13.73+1,14.28*fs]);

DS = 1; % down sample further if requested
D = 40; % number channels (don't set too high)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
soundPath = '../../tsgp/datasets/audio/';
[y,fs] = wavread([soundPath,File,'.wav']); % reads in the file
yTrain = y(RngLimTrain(1):RngLimTrain(2),1); 
yTrain = resample(yTrain, fs/DS, fs); % downsample the input
fs = fs/DS;
yTrain = yTrain/sqrt(var(yTrain)); % rescale the input to unit variance
T = length(yTrain);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Setup the filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
trainFilter = 0;

%if trainFilter==1
  % Learn properties of the filters (centre frequency and width)
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 500;
  opts.minT = 5000;
  opts.numIts = 10;
  DInit = D;
  [Var1Fit,Lam1Fit,omFit,InfoFit] = fit_probSTFT(yTrain,DInit,opts); % trains filters to
                                              % match the spectrum

  % Order carriers by centre frequency
  [omFit,ind] = sort(omFit);
  Var1Fit = Var1Fit(ind); 
  Lam1Fit = Lam1Fit(ind); 

  % useful to know the bandwidths and marginal variances of the
  % carriers, so computing them here:

  [fmax,df,varMa] = probSpec2freq(omFit,Lam1Fit,Var1Fit);

  %ySampNoise = samplePFB(Lam1,Var1,om,0,T);
%else
  dfFrac = 1/10; % best results dfFrac = 1/20
  fmax = logspace(log10(1/100),log10(0.2),D)';
  [om,Lam1,Var1] = freq2probSpec(fmax,fmax*dfFrac,ones(D,1)/D);
%end

% plots the sprectrum of the filter
figH1 = plot_pSTFT(Var1Fit,omFit,Lam1Fit,fs,1);
figH2 = plot_pSTFT(Var1,om,Lam1,fs,1);