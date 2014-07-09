clear;

% runs NMF on a speech sound

% Specify where to load the data from
soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';

% Specify where to save the data to
saveDir = '~/data/probFB/nmf/';


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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply the filter and produce the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vary = 0;
ZTrain = probFB(yTrain,Lam1,Var1,om,0); % applies filters to the signal, replaced AR2 filter bank (much faster)

ATrain = abs(ZTrain').^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = 40;%D;  % number of features

HInit = exp(randn(T,K));
ks = ceil(T*rand(K,1));
WInit = ATrain(ks,:);
vary = zeros(T,D);

% initialise with regular NMF
Opts.restarts = 50;
Opts.numIts = 5000;
%[WEst1,HEst1,info1] = nmf(ATrain,WInit,HInit,[],[],[],vary,Opts);
[WEst1,HEst1,info1] = nmf_fp(ATrain,WInit,HInit,vary,Opts);

% order according to slowness (mean square derivative)
fastness = mean(diff(HEst1).^2)./var(HEst1);
[val,ind] = sort(fastness,'descend');
HEst1 = HEst1(:,ind); 
WEst1 = WEst1(ind,:);


% % TRUNCATED GAUSSIAN PRIORS
% % trying truncated Gaussian temporal priors instead

% % set lam from the correlation between successive hs
% muinf = mean(HEst1);
% varinf = var(HEst1);
% lam = mean((HEst1(1:T-1,:)-ones(T-1,1)*muinf).*(HEst1(2:T,:)- ones(T-1,1)*muinf))./var(HEst1);

% Opts.numIts = 20000;
% Opts.restarts = 1;
% [WEst2,HEst2,info2] = nmf(ATrain,WEst1,HEst1,[],varinf,lam,vary,Opts);
% %[WEst2,HEst2,info2] = nmf(ATrain,WEst2,HEst2,[],varinf,lam,vary,Opts);

%% LOG-NORMAL SE GAUSSIAN PROCESS TEMPORAL NMF
% initialise

lenx = zeros(K,1);
mux = zeros(K,1);
varx = zeros(K,1);

for k=1:K
  % threshold
  logHthresh = log(HEst1(:,k)+1e-8);
  
  % filter
  filt = exp(-1/2*([-100:100].^2)/(1^2));
  filt = filt/sum(filt);
  logHsm = conv(logHthresh,filt,'same');
  
  % fit GP
  mux(k) = mean(logHsm);
  [lenx(k),varx(k),info] = trainSEGP_RS(logHsm-mux(k));
  
%  figure
%  hold on
%  plot(logHthresh,'-b')
%  plot(logHsm,'-k')
%  keyboard
end

disp(['length scale parameters'])
lenx

Opts.numIts = 200;
[HEst2,info2] = tnmf_inf(ATrain,WEst1,HEst1+1e-8,lenx,mux,varx,vary,Opts);

Opts.numIts = 2000;
[WEst2,HEst2,info2] = tnmf(ATrain,WEst1,HEst2,lenx,mux,varx,vary,Opts);

OptsJoint.numIts = 700; % 50 iterations will take 20mins
OptsJoint.progress_chunk = 100;
varyGTF = 0.01;
HEst2Mod = HEst2+1e-5;
WEst2Mod = WEst2+1e-5;

lenx2 = 100*lenx;% I think the length scales should be longer than: 1./(fmax*dfFrac)
mux2 = mux;
varx2 = varx;

% For initialising using broad features
% sigf = 3;
% ind = ones(D,1)*[1:D];
% dind = (ind-ind').^2;
% WEst2Mod = exp(-1/(2*sigf^2)*(dind.^2))+1e-5;
% lenx2 = 100*mean(lenx)*ones(D,1);
% mux2 = mean(mux)*ones(D,1);
% varx2 = mean(varx)*ones(D,1);
% HEst2Mod = exp(randn(T,D)-6);

[WEst3,HEst3,mnV3,covV3,info3] = GTFtNMF(yTrain,WEst2Mod,HEst2Mod,Lam1,Var1,om,varyGTF, ...
				       lenx2,mux2,varx2,OptsJoint);

clear mnV3 covV3
%[WEst3,HEst3,mnV3,covV3,info3] = GTFtNMF(yTrain,WEst3,HEst3,Lam1,Var1,om,varyGTF,lenx2,mux2,varx2,OptsJoint);
%mux2 = mux2+6;

%save('/home/ret26/data/probFB/nmf/trainDenoise9_wide_D40_med.mat')

filename = 'trainDenoise10_D40';
load(['/home/ret26/data/probFB/nmf/',filename,'.mat'])

OptsJointSpec.numIts = 100; % 50 iterations will take 20mins
OptsJointSpec.progress_chunk = 25;

[WEst4,HEst4,Lam14,Var14,om4,mnV4,covV4,info4] = GTFtNMF_train_spec(yTrain, ...
						  WEst3,HEst3,Lam1, ...
						  Var1,om,varyGTF, ...
						  lenx2,mux2,varx2,OptsJointSpec);

[WEst4,HEst4,Lam14,Var14,om4,mnV4,covV4,info4] = GTFtNMF_train_spec(yTrain, ...
						  WEst4,HEst4,Lam14, ...
						  Var14,om4,varyGTF, ...
						  lenx2,mux2,varx2,OptsJointSpec);

clear mnV4 covV4

filename = 'trainDenoise11_D40';
save(['/home/ret26/data/probFB/nmf/',filename,'.mat'])

figH2 = plot_pSTFT(Var1,om,Lam1,fs,0);
figH3 = plot_pSTFT(Var14,om4,Lam14,fs,0);
filename = 'trainDenoise10_D40';


% PLOT RESULTS
figure
hold on
plot(info1.Obj,'-k')
plot(info2.Obj,'-r')
plot(info3.Obj,'-m')
 

% for k=1:K
%   %subplot(K,1,k)
%   figure
%   hold on
%   plot(HEst1(:,k),'-k')
%   plot(HEst2(:,k),'-r')
%   plot(HEst3(:,k),'-m')
  
%   set(gca,'yscale','log')
%   set(gca,'ylim',[min(HEst3(:,k)),max(HEst3(:,k))+1])
% end

figure
subplot(3,1,1)
imagesc(log(ATrain)')

subplot(3,1,2)
imagesc(log(HEst1*WEst1)')

subplot(3,1,3)
imagesc(log(HEst2*WEst2)')

%mean(abs(HEst1*WEst1-ATrain))
%mean(abs(HEst2*WEst2-ATrain))
snr1_a = 10*(log10(mean(ATrain.^2,1))-log10(mean((HEst1*WEst1-ATrain).^2,1)));
snr2_a = 10*(log10(mean(ATrain.^2,1))-log10(mean((HEst2*WEst2-ATrain).^2,1)));

snr1_loga = 10*(log10(mean(log(ATrain).^2,1))-log10(mean((log(HEst1*WEst1)-log(ATrain)).^2,1)));
snr2_loga = 10*(log10(mean(log(ATrain).^2,1))-log10(mean((log(HEst2*WEst2)-log(ATrain)).^2,1)));

disp(['upper limit on SNRs ',num2str([mean(snr2_a),mean(snr2_loga)])]);

figure
subplot(2,2,1)
hold on
title('NMF')
imagesc(WEst1)
set(gca,'xlim',[1,D],'ylim',[1,K])

subplot(2,2,2)
hold on
title('tNMF')
imagesc(WEst2)
set(gca,'xlim',[1,D],'ylim',[1,K])

subplot(2,2,3)
hold on
title('NMF')
imagesc(log(WEst1))
set(gca,'xlim',[1,D],'ylim',[1,K])

subplot(2,2,4)
hold on
title('tNMF')
imagesc(log(WEst2))
set(gca,'xlim',[1,D],'ylim',[1,K])


figure
subplot(2,2,1)
hold on
title('tNMF')
imagesc(WEst2)
set(gca,'xlim',[1,D],'ylim',[1,K])

subplot(2,2,2)
hold on
title('GTFtNMF')
imagesc(WEst3)
set(gca,'xlim',[1,D],'ylim',[1,K])

subplot(2,2,3)
hold on
title('tNMF')
imagesc(log(WEst2))
set(gca,'xlim',[1,D],'ylim',[1,K])

subplot(2,2,4)
hold on
title('GTFtNMF')
imagesc(log(WEst3))
set(gca,'xlim',[1,D],'ylim',[1,K])
