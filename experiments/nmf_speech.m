clear;

% runs NMF on a speech sound

% Specify where to load the data from
soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';

% Specify where to save the data to
saveDir = '~/Data/probFB/nmf/';

% load signal
File = '74 - Sentences'; % Name of file to load
fs = 16000; % sampling rate of file
RngLimTrain = round([fs*1/2+1,2.1*fs]);  % Picks a small range
RngLimTest = round([fs*13.7+1,15.4*fs]);  % Picks a small range

DS = 1; % down sample further if requested
D = 20; % number channels (don't set too high)

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

if trainFilter==1
  % Learn properties of the filters (centre frequency and width)
  opts.verbose = 1; % view plots of the fitting process
  opts.minT = 500;
  opts.minT = 5000;
  opts.numIts = 10;
  DInit = D;
  [Var1,Lam1,om,Info] = fit_probSTFT(yTrain,DInit,opts); % trains filters to
                                              % match the spectrum

  % Order carriers by centre frequency
  [om,ind] = sort(om);
  Var1 = Var1(ind); 
  Lam1 = Lam1(ind); 

  % useful to know the bandwidths and marginal variances of the
  % carriers, so computing them here:

  [fmax,df,varMa] = probSpec2freq(om,Lam1,Var1);

  %ySampNoise = samplePFB(Lam1,Var1,om,0,T);
else
  dfFrac = 1/15;
  fmax = logspace(log10(1/50),log10(0.3),D)';
  [om,Lam1,Var1] = freq2probSpec(fmax,fmax*dfFrac,ones(D,1)/D);
end

% plots the sprectrum of the filter
figH = plot_pSTFT(Var1,om,Lam1,fs,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Apply the filter and produce the spectrogram
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vary = 0;
ZTrain = probFB(yTrain,Lam1,Var1,om,0); % applies filters to the signal, replaced AR2 filter bank (much faster)

ATrain = abs(ZTrain').^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Train the model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
K = D;  % number of features

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

% PLOT RESULTS
figure
hold on
plot(info1.Obj,'-r')
plot(info2.Obj,'-k')
 
figure
for k=1:K
  subplot(K,1,k)
  hold on
  plot(HEst1(:,k),'-r')
  plot(HEst2(:,k),'-k')

  
set(gca,'yscale','log')
end

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

% subplot(2,1,1)
% imagesc(W)
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TESTING


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Load signal and pre-process
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[y,fs] = wavread([soundPath,File,'.wav']); % reads in the file
yTest = y(RngLimTest(1):RngLimTest(2),1); 
yTest = resample(yTest, fs/DS, fs); % downsample the input
fs = fs/DS;
yTest = yTest/sqrt(var(yTest)); % rescale the input to unit variance
T = length(yTest);

% clean spectrogram computed here
ZTest = probFB(yTest,Lam1,Var1,om,0); 
ATest = abs(ZTest').^2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Denoising
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

L = 3;
varys = [logspace(log10(1e-1),log10(40),L)];

%L=1;
%varys = 1;

HInit = exp(randn(T,K))/10000;

mse_a = zeros(L,1);
mse_loga = zeros(L,1);
mse_orig_a = zeros(L,1);
mse_orig_loga = zeros(L,1);
snr_orig_a = zeros(L,D);
snr_a = zeros(L,D);
snr_orig_loga = zeros(L,D);
snr_loga = zeros(L,D);
snr_orig_y = zeros(L,1);
snr_y = zeros(L,1);
pesq_before = zeros(L,1);
pesq_after = zeros(L,1);


opts_inf.numIts = 1000;
opts_inf.progress_chunk = 100;
%opts_inf.restarts = 0;

for l=1:L
  disp(['Progress ',num2str(l),'/',num2str(L)]);

  % set the state 
  randn('state',1);
  
  % noisy signal
  yNoisy = yTest + randn(T,1)*sqrt(varys(l));
  ZNoisy = probFB(yNoisy,Lam1,Var1,om,0); 
  ANoisy = abs(ZNoisy').^2;

  % figure out the noise levels -- this is done empirically, but we
  % could do this analytically if needed
  filt_noise = probFB(randn(T,1)*sqrt(varys(l)),Lam1,Var1,om,0);
  varyCur = ones(T,1)*var(filt_noise');
  
  % denoise using NMF
%  [HTest,info] =
%  nmf_inf(ANoisy,WEst2,HInit,[],[],[],varyCur,opts_inf);


%  [HTest_fp,infoInit] = nmf_inf_fp(ANoisy,WEst2,HInit,varyCur);
%  [HTest,info] = tnmf_inf(ANoisy,WEst2,HTest_fp+1e-9,lenx,mux,varx,varyCur,opts_inf);

[HTest,info] = tnmf_inf(ANoisy,WEst2,HInit,lenx,mux,varx,varyCur,opts_inf);

% truncated Gaussian 
%  [HTest_fp,infoInit] = nmf_inf_fp(ANoisy,WEst2,HInit,varyCur);
%  [HTest,info] = nmf_inf(ANoisy,WEst2,HTest_fp,[],varinf,lam,varyCur,opts_inf);


  ADenoise = HTest*WEst2;
  
  % figure out the amount NMF has denoised the spectrogram by
  mse_orig_a(l) = mean((ATest(:)-ANoisy(:)).^2);
  mse_orig_loga(l) = mean((log(ATest(:))-log(ANoisy(:))).^2);
  mse_a(l) = mean((ATest(:)-ADenoise(:)).^2);
  mse_loga(l) = mean((log(ATest(:))-log(ADenoise(:))).^2);

  snr_orig_a(l,:) = 10*(log10(sum(ATest.^2))-log10(sum((ATest-ANoisy).^2)));
  snr_orig_loga(l,:) = 10*(log10(sum(log(ATest).^2))-log10(sum((log(ATest)-log(ANoisy)).^2)));
  snr_a(l,:) = 10*(log10(sum(ATest.^2))-log10(sum((ATest-ADenoise).^2)));
  snr_loga(l,:) = 10*(log10(sum(log(ATest).^2))-log10(sum((log(ATest)-log(ADenoise)).^2)));

  % reconstruction of the signal via least squares
  spec = get_probFB_spec(Lam1,Var1,om,0,T);
  lamRec = 0; varyRec = (1-lamRec^2);
  yInit = randn(T,1)/1000;
  [yDenoise,aRec,info3] = recon_FB_mag(yInit,ADenoise.^0.5,spec,[],lamRec,varyRec);

  % evaluate SNR and perceptual improvement
  snr_orig_y(l) = 10*(log10(sum(yTest.^2))-log10(sum((yTest-yNoisy).^2)));
  snr_y(l) = 10*(log10(sum(yTest.^2))-log10(sum((yTest-yDenoise).^2)));

  pesq_before(l) = pesq(yTest, yNoisy, fs);
  pesq_after(l) = pesq(yTest, yDenoise, fs);

  
  figure
  subplot(3,1,1)
  imagesc(log(ATest)')
  
  subplot(3,1,2)
  imagesc(log(ADenoise)')
  
  subplot(3,1,3)
  imagesc(log(ANoisy)')
%  keyboard
end

figure
subplot(1,2,1)
hold on
title('spectrogram')
plot(mean(snr_orig_a,2),mean(snr_a-snr_orig_a,2),'.-')
xlabel('mean SNR before /dB')
ylabel('mean SNR improvement /dB')

subplot(1,2,2)
hold on
title('log-spectrogram')
plot(mean(snr_orig_loga,2),mean(snr_loga-snr_orig_loga,2),'.-')
xlabel('mean SNR before /dB')
ylabel('mean SNR improvement /dB')


figure
subplot(2,1,1)
hold on
title('waveform')
plot(pesq_before,pesq_after-pesq_before,'.-')
xlabel('PSEQ before')
ylabel('PSEQ improvement')

subplot(2,1,2)
hold on
plot(snr_orig_y,snr_y-snr_orig_y,'.-')
xlabel('SNR before /dB')
ylabel('SNR improvement /dB')

[mean(snr_a-snr_orig_a,2)';
mean(snr_loga-snr_orig_loga,2)']
[pesq_before,pesq_after-pesq_before]
[snr_orig_y,snr_y-snr_orig_y]


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %% GAP RESTORATION
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% L = 3;
% gapLim = [5,100];
% NGaps = 20;
% gaps = ceil(logspace(log10(gapLim(1)),log10(gapLim(2)),L));
% gapPos = gapLim(2)+ceil(rand(NGaps,1)*(T-2*gapLim(2)));
  
% HInit = exp(randn(T,K))/10000;

% recon_mse_a = zeros(L,1);
% recon_mse_loga = zeros(L,1);
% recon_mse_orig_a = zeros(L,1);
% recon_mse_orig_loga = zeros(L,1);
% recon_snr_orig_a = zeros(L,D);
% recon_snr_a = zeros(L,D);
% recon_snr_orig_loga = zeros(L,D);
% recon_snr_loga = zeros(L,D);
% recon_snr_orig_y = zeros(L,1);
% recon_snr_y = zeros(L,1);
% recon_pesq_before = zeros(L,1);
% recon_pesq_after = zeros(L,1);

% opts_inf.numIts = 1000;
% opts_inf.progress_chunk = 100;

% for l=1:L
%   yGap = yTest;
%   varyCur = zeros(T,D);
%   ind = [];
%   for ng=1:NGaps
%     ind = [ind,gapPos(ng)+[-ceil(gaps(l)/2):+ceil(gaps(l)/2)]];
%   end
  
%   yGap(ind) = 0;
%   varyCur(ind,:) = 1e5;
  
%   ZGap = probFB(yGap,Lam1,Var1,om,0); 
%   AGap = abs(ZGap').^2;

%   % inference
%   [HRecon,info] = tnmf_inf(AGap,WEst2,HInit,lenx,mux,varx,varyCur,opts_inf);

%   ARecon = HRecon*WEst2;
  
%   % figure out the amount NMF has denoised the spectrogram by
%   recon_mse_orig_a(l) = mean(mean((ATest(ind,:)-AGap(ind,:)).^2));
%   recon_mse_orig_loga(l) = mean(mean((log(ATest(ind,:))-log(AGap(ind,:))).^2));
%   recon_mse_a(l) = mean(mean((ATest(ind,:)-ARecon(ind,:)).^2));
%   recon_mse_loga(l) = mean(mean((log(ATest(ind,:))-log(ARecon(ind,:))).^2));

%   recon_snr_orig_a(l,:) = 10*(log10(sum(ATest(ind,:).^2))-log10(sum((ATest(ind,:)-AGap(ind,:)).^2)));
%   recon_snr_orig_loga(l,:) = 10*(log10(sum(log(ATest(ind,:)).^2))-log10(sum((log(ATest(ind,:))-log(AGap(ind,:))).^2)));
%   recon_snr_a(l,:) = 10*(log10(sum(ATest(ind,:).^2))-log10(sum((ATest(ind,:)-ARecon(ind,:)).^2)));
%   recon_snr_loga(l,:) = 10*(log10(sum(log(ATest(ind,:)).^2))-log10(sum((log(ATest(ind,:))-log(ARecon(ind,:))).^2)));

%   % reconstruction of the signal via least squares
%   spec = get_probFB_spec(Lam1,Var1,om,0,T);
%   lamRec = 0; varyRec = (1-lamRec^2);
%   yInit = randn(T,1)/1000;
%   [yRecon,aRec,info3] = recon_FB_mag(yInit,ARecon.^0.5,spec,[],lamRec,varyRec);

%   % evaluate SNR and perceptual improvement
%   recon_snr_orig_y(l) = 10*(log10(sum(yTest(ind).^2))-log10(sum((yTest(ind)-yGap(ind)).^2)));
%   recon_snr_y(l) = 10*(log10(sum(yTest(ind).^2))-log10(sum((yTest(ind)-yRecon(ind)).^2)));

%   recon_pesq_before(l) = pesq(yTest, yGap, fs);
%   recon_pesq_after(l) = pesq(yTest, yRecon, fs);

  
%   figure
%   subplot(3,1,1)
%   imagesc(log(ATest)')
  
%   subplot(3,1,2)
%   imagesc(log(ARecon)')
  
%   subplot(3,1,3)
%   imagesc(log(AGap)')
%   keyboard

% end



% figure
% subplot(1,2,1)
% hold on
% title('spectrogram')
% plot(mean(recon_snr_orig_a,2),mean(recon_snr_a-recon_snr_orig_a,2),'.-')
% xlabel('mean SNR before /dB')
% ylabel('mean SNR improvement /dB')

% subplot(1,2,2)
% hold on
% title('log-spectrogram')
% plot(mean(recon_snr_orig_loga,2),mean(recon_snr_loga-recon_snr_orig_loga,2),'.-')
% xlabel('mean SNR before /dB')
% ylabel('mean SNR improvement /dB')


% figure
% subplot(2,1,1)
% hold on
% title('waveform')
% plot(recon_pesq_before,recon_pesq_after-recon_pesq_before,'.-')
% xlabel('PSEQ before')
% ylabel('PSEQ improvement')

% subplot(2,1,2)
% hold on
% plot(recon_snr_orig_y,recon_snr_y-recon_snr_orig_y,'.-')
% xlabel('SNR before /dB')
% ylabel('SNR improvement /dB')

% [mean(recon_snr_a-recon_snr_orig_a,2)';
% mean(recon_snr_loga-recon_snr_orig_loga,2)']
% [recon_pesq_before,recon_pesq_after-pesq_before]
% [recon_snr_orig_y,recon_snr_y-recon_snr_orig_y]

