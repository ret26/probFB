function test_suite = test_recon_FB_mag
  initTestSuite;

% Tests:
%
% function [y,info] = recon_FB_mag(y,aTar,spec,varargin)
% 
  
function test_reconstruct_speech

% shows that the function does sensible things for a speech example
% for this example the best performing reconstruction method is to use
% a little bit of regularisation (although with no temporal
% constraint) and to use the squared error in the envelopes as the
% error-metric. For this example, the pesq error introduced by this
% process is about equivalent to adding Gaussian noise of variance = 0.05

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOAD DATA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Specify where to load the data from
soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';

% Specify where to save the data to
saveDir = '~/Data/probFB/nmf/';

% load signal
File = '74 - Sentences'; % Name of file to load
fs = 16000; % sampling rate of file
RngLim = round([fs*1/2+1,2.1*fs]);  % Picks a small range
%RngLim = round([fs*1/2+1,1*fs]);  % Picks a small range
DS = 1; % down sample further if requested
D = 5; % number channels (don't set too high)

[y,fs] = wavread([soundPath,File,'.wav']); % reads in the file
y = y(RngLim(1):RngLim(2),1); 
y = resample(y, fs/DS, fs); % downsample the input
fs = fs/DS;
y = y/sqrt(var(y)); % rescale the input to unit variance
T = length(y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% APPLY FILTERBANK -- pSPEC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fmax = logspace(log10(0.01),log10(0.4),D);
df = fmax/10;
varMa = ones(1,D);
[om,lamx,varx] = freq2probSpec(fmax,df,varMa);

vary = 0;
Z = probFB(y,lamx,varx,om,vary); % applies filters to the signal, replaced AR2 filter bank (much faster)

aTar = abs(Z)';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% RECONSTRUCT 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

spec = get_probFB_spec(lamx,varx,om,vary,T);

yInit = randn(T,1)/1000;

% vanilla squared error in envelopes, no regularisation
[yRec1,aRec1,info1] = recon_FB_mag(yInit,aTar,spec);

% log-error
% (didn't perform as well as differences in the envelopes)
opts.error_measure = 'loga';
[yRec2,aRec2,info2] = recon_FB_mag(yInit,aTar,spec,opts);

% regularisation 1
% no temporal regularisation
opts.error_measure = 'a';
opts.numIts =1000; % higher numbers of iterations don't make much
                   % difference to the perceptual quality rating
opts.progress_chunk = 200;
lam = 0;
vary = (1-lam^2);
[yRec3,aRec3,info3] = recon_FB_mag(yInit,aTar,spec,opts,lam,vary);
%[yRec3,aRec3,info3] = recon_FB_mag(yInit,aTar,spec,[],lam,vary);

% regularisation 2
lam = 0.1;
vary = 1-lam^2;
yInit = randn(T,1)/1000;
[yRec4,aRec4,info4] = recon_FB_mag(yInit,aTar,spec,[],lam,vary);

% evaluate perceptual score of reconstructions
ratB = pesq(y, yInit, fs);
rat1 = pesq(y, yRec1, fs);
rat2 = pesq(y, yRec2, fs);
rat3 = pesq(y, yRec3, fs);
rat4 = pesq(y, yRec4, fs);


disp(['----------------------'])
disp(['before ',num2str(ratB)]);
disp(['vanilla ',num2str(rat1)]);
disp(['log-error ',num2str(rat2)]);
disp(['reg lam=0 ',num2str(rat3)]);
disp(['reg lam=0.1 ',num2str(rat4)]);
disp(['----------------------'])

% compare the reconstruction perceptual scores to corrupting the
% signal with Gaussian noise
L = 40;
varys = linspace(0,1,L);
rat = zeros(L,1);
for l=1:L
  rat(l) = pesq(y, y+randn(T,1)*sqrt(varys(l)), fs);
end

[valB,locB] = min(abs(ratB-rat));
[val1,loc1] = min(abs(rat1-rat));
[val2,loc2] = min(abs(rat2-rat));
[val3,loc3] = min(abs(rat3-rat));
[val4,loc4] = min(abs(rat4-rat));

figure
hold on
plot(varys,rat,'-k')
plot(varys(locB),ratB,'og')
plot(varys(loc1),rat1,'or')
plot(varys(loc2),rat2,'oc')
plot(varys(loc3),rat3,'ok')
plot(varys(loc4),rat4,'ob')
legend('noise corrupted signal','before','no reg SE a','no reg SE loga','reg lam=0 SE a','reg lam=0.1 SE a')


figure
hold on;
plot(y,'-k')
plot(yRec3,'-r')

figure
subplot(2,1,1)
imagesc(log(aTar'))

subplot(2,1,2)
imagesc(log(aRec3'))

keyboard

%tol = 1e-6;
%assertElementsAlmostEqual(vals,0,'absolute',tol,0)
