clear;

% runs NMF on a speech sound

% Specify where to load the data from
soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';

% Specify where to save the data to
saveDir = '~/Data/probFB/nmf/';

% load signal
File = '74 - Sentences'; % Name of file to load
fs = 16000; % sampling rate of file
RngLim = round([fs*1/2+1,2.1*fs]);  % Picks a small range
DS = 1; % down sample further if requested
D = 10; % number channels (don't set too high)
K = D;  % number of features

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,fs] = wavread([soundPath,File,'.wav']); % reads in the file
y = y(RngLim(1):RngLim(2),1); 
y = resample(y, fs/DS, fs); % downsample the input
fs = fs/DS;
y = y/sqrt(var(y)); % rescale the input to unit variance
T = length(y);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Train filter
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Learn properties of the filters (centre frequency and width)
opts.verbose = 1; % view plots of the fitting process
[Var1,Lam1,om,Info] = fit_probSTFT(y,D,opts); % trains filters to
                                              % match the spectrum
% Order carriers by centre frequency
[om,ind] = sort(om);
Var1 = Var1(ind); 
Lam1 = Lam1(ind); 

% useful to know the bandwidths and marginal variances of the
% carriers, so computing them here:

[fmax,df,varMa] = probSpec2freq(om,Lam1,Var1);

ySampNoise = samplePFB(Lam1,Var1,om,0,T);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Filtering step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vary = 0;
Z = probFB(y,Lam1,Var1,om,vary); % applies filters to the signal, replaced AR2 filter bank (much faster)

A = abs(Z)';

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Forward model Settings

% muinf = mean(A(:));
% varinf = var(A(:));
% lam = linspace(0.1,0.9,K);
% vary = zeros(T,1);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % TEST

% WInit = exp(randn(K,D))/1;
% HInit = exp(randn(T,K))/1;

% Opts.restarts = 5;
% [WEst,HEst,info] = nmf(A,WInit,HInit,muinf,varinf,lam,vary,Opts);

% figure
% plot(info.Obj)

% figure
% for k=1:K
%   subplot(K,1,k)
%   hold on
%   plot(H(:,k),'-k')
%   plot(HEst(:,k),'-r')
% set(gca,'yscale','log')
% end

% figure
% subplot(2,1,2)
% imagesc(WEst)

% subplot(2,1,1)
% imagesc(W)

