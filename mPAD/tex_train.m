clear;

% add all of the necessary paths
path1 = genpath('~/projects/code/imageModels');
path2 = genpath('~/projects/code/pad');
path3 = genpath('~/projects/code/sfa');
path4 = genpath('~/projects/code/probFB');
path5 = genpath('~/projects/code/conj-grad');
path6 = genpath('~/projects/code/matlab_xunit');

addpath(path1);
addpath(path2);
addpath(path3);
addpath(path4);
addpath(path5);
addpath(path6);

% Specify where to load the data from
soundPath = '~/Music/Various Artists/Best Ever Sound Effects - Vol.3 - Sounds Of Nature Sound Effects/';

% Specify where to save the data to
saveDir = '~/sandbox/kyriacos/sounds/';


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Select signal
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Explanation of variables:
% File =  Name of file to load
% fs =  sampling rate of file
% RngLim =  Range of file to train off - usually about 5s is long enough
% DS = down sample by this factor (DS=1 => no down sampling)
% D =  number channels (don't set too high)
% K = number of features to retain (could look at the eigValsPCA to decide how many components to retain)
% fac = initialise the time-scales of the modulation using fac/channel-centre-frequency


% File = '83 - Stream'; % Name of file to load
% fs = 48000;
% RngLim = round([fs*1/2+1,5.1*fs]);  % Picks a small range
% DS = 2; % should drop to DS = 2 eventually
% fac = 4;
% D = 20; % number channels (don't set too high)
% K = D;

% File = '24 - Fire Burning'; % Name of file to load
% fs = 44100;
% RngLim = round([fs*3+1,7.1*fs]);  % Picks a small range
% DS = 2; % should drop to DS = 2 eventually
% fac = 4;
% D = 20; % number channels (don't set too high)
% K = D;

% File = '41 - Mountain Stream - Deep'; % Name of file to load
% fs = 44100;
% RngLim = round([fs*3+1,10.1*fs]);  % Picks a small range
% DS = 2; % should drop to DS = 2 eventually
% fac = 4;
% D = 20; % number channels (don't set too high)
% K = D;

% File = '05 - Bubbling stream'; % Name of file to load
% fs = 44100;
% %RngLim = round([fs*3+1,8.1*fs]);  % Picks a small range
% RngLim = round([fs*3+1,8.1*fs]);  % Picks a small range
% DS = 2; % should drop to DS = 2 eventually
% fac = 4;
% D = 20; % number channels (don't set too high)
% K = D;

% File = '80 - Fire'; % Name of file to load
% fs = 10989;
% RngLim = round([fs*0+1,3.5*fs]);  % Picks a small range
% DS = 1; % should drop to DS = 2 eventually
% fac = 4;
% D = 20; % number channels (don't set too high)
% K = D;

% File = '85 - Applause'; % Name of file to load
% fs = 44100;
% RngLim = round([fs*10+1,floor(17*fs)]);  % Picks a small range
% DS = 2; % should drop to DS = 2 eventually
% fac = 2;
% D = 20; % number channels (don't set too high)
% K = 3;

%%%%%%%%%%%%
% Training for the following sounds is not working very well at present:

% File = '75 - Wind'; % Name of file to load
% fs = 11025;
% RngLim = round([fs*0+1,floor(3.9*fs)]);  % Picks a small range
% DS = 1; % should drop to DS = 2 eventually
% fac = -1500;
% D = 20; % number channels (don't set too high)
% K = 3;

% File = '81 - Wind'; % Name of file to load
% fs = 11025;
% RngLim = round([fs*0+1,floor(7.5*fs)]);  % Picks a small range
% DS = 1; % should drop to DS = 2 eventually
% opts.minT = 1000; opts.maxT = 2000;
% fac = -1000;
% D = 20; % number channels (don't set too high)
% K = 3;

% File = '54 - Rain On Leaves'; % Name of file to load
% fs = 44100;
% %RngLim = round([fs*13.5+1,floor(18.5*fs)]);  % Picks a small range
% RngLim = round([fs*1+1,floor(11*fs)]);  % Picks a small range
% DS = 2; 
% fac = 2;
% D = 20; % number channels (don't set too high)
% K = 15;

File = '74 - Sentences'; % Name of file to load
fs = 16000; % sampling rate of file
RngLim = round([fs*1/2+1,2.1*fs]);  % Picks a small range
DS = 1; % down sample further if requested
D = 20; % number channels (don't set too high)
K = D;  % number of features
fac = 4; % initialise the time-scales of the modulation using fac/channel-centre-frequency

% File = '88 - Rain Noisy';
% fs = 44100;
% RngLim = round([fs*3+1,floor(6*fs)]);  % Picks a small range
% DS = 2; 
% fac = 2;
% D = 20; % number channels (don't set too high)
% K = 10;

% File = '89 - Fire'; % Name of file to load
% fs = 44100;
% RngLim = round([fs*3+1,10.1*fs]);  % Picks a small range
% DS = 2; % should drop to DS = 2 eventually
% fac = 3;
% D = 20; % number channels (don't set too high)
% K = 3;

% File = '90 - Fire'; % Name of file to load
% fs = 22050;
% RngLim = round([fs*3+1,10.1*fs]);  % Picks a small range
% DS = 1; 
% fac = 4;
% D = 20; % number channels (don't set too high)
% K = 3;

% File = '91 - Rain'; % Name of file to load
% fs = 44100;
% %RngLim = round([fs*13.5+1,floor(18.5*fs)]);  % Picks a small range
% RngLim = round([fs*1+1,floor(8*fs)]);  % Picks a small range
% DS = 2; 
% fac = 2;
% D = 20; % number channels (don't set too high)
% K = 10;

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

keyboard

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Demodulation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% time-scale in samples that we demodulate at (depends on the
% centre frequency - higher frequency filters support faster modulation)

if fac>0
  len = fac*(2*pi./om)';
else
  len = -fac*ones(1,D);
end


% demodulate using probabilistic amplitude demodulation
[A,C,ParamsGPPAD,InfoGPPAD] = GPPAD(real(Z'),len,1); 

%% STFT code if needed:
% S = zeros(D,T);
% for d=1:D
%   S(d,:) = exp(-i*om(d)*[1:T]).*S(d,:);
% end

%% there's some sorting out of scales to do here
A = A'./(sqrt(ParamsGPPAD.varc)*ones(1,T)); % rescaling A
%A = A'.*/(sqrt(var(C,2))*ones(1,T));

C = real(Z)./A; % recomputing the carriers so that X1=A.*C

d=1; figure; hold on; plot(real(Z(d,:)),'-k'); plot(A(d,:),'-r','linewidth',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Pre-process & reduce dimensionality before feature extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


X2 = log(exp(A')-1); % convert positive envelopes into +/- quantities

muX = mean(X2); 
X2 = X2 - ones(T,1)*muX; % remove the mean

[WPCA,XPCA,eigValsPCA] = PCA(X2); % PCA for dimensionality reduction
GPCA = inv(WPCA);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Feature extraction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% One option is to use SFA

[WSFA,XSFA,omSqSFA] = SFA(XPCA(:,1:K));
GSFA = inv(WSFA);
G1 = GSFA*GPCA(1:K,:); %log(A') \approx XSFA*G1+ones(T,1)*muX;

%Note that: omSqSFA = (2*pi)^2*getVarSpec(XSFA)
LenSFA = 1./sqrt(omSqSFA)'; % time-scales of the SFA
                            % solution in samples

VarSFA = var(XSFA); % variance of the SFA solutions (will be = 1)

% Another option is to use ICA 
% this will probably result in shorter time-scales as it is not
% explicitly maximising slowness like SFA

NumIts = 2000; % set the number of iterations. Choose a value for which
               % the objective has converged - check convergence via
               % plotObj(infoICA.Obj)
	       
[WICA,XICA,infoICA] = studentICA(XPCA(:,1:K),NumIts);

GICA = inv(WICA)'; %NB the transpose here because the ICA weight
                   %matrix is defined the other way around

G2 = GICA*GPCA(1:K,:); % log(A') \approx XICA*G2+ones(T,1)*muX;

omSqICA = (2*pi)^2*getVarSpec(XICA);
LenICA = 1./sqrt(omSqICA)'; % time-scales of the ICA
			   % solution in samples
VarICA = var(XICA); % variance of the ICA solutions

% NB. To view the features extracted by the above methods in
% log-envelope space: use plot(G2(k,:)) to see the spectral shape or
% imagesc(XICA(:,k)*G2(k,:)+ones(T,1)*muX) to see the activity pattern
% resulting from that feature.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Fine-tune
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Either initialise using ICA or SFA weights
ICA = 0; % 0 => use SFA, 1 => use ICA (SFA seems a little better at
         % present so I've used that by default)

% Carrier properties
Var1 = var(C')'.*(1-Lam1.^2);%need to rescale the carrier energy
Mu2 = muX';
vary = 0.001;

% ICA
ParamsICA = packParamsMPAD(G2',Lam1,Var1,LenICA,VarICA,Mu2,vary,om); 

% SFA
ParamsSFA = packParamsMPAD(G1',Lam1,Var1,LenSFA,VarSFA,Mu2,vary,om); 

if ICA==1
  Params = ParamsICA;
  X2 = XICA;
  typ = 'ICA';
else
  Params = ParamsSFA;
  Params.G = ParamsSFA.G*1.4; % ParamsSFA.G*1.4 - seems to make things sound much
                                 % better presumably countering
                                 % biases in the amplitude estimation
  X2 = XSFA;
  typ = 'SFA';
  figure
  imagesc(ParamsSFA.G)
end

optsCG.numIts = 50;
optsCG.progress_chunk = 5;

opts.numIts = 200;
opts.progress_chunk = 5;
%[X1FT,X2FT,ParamsFT] = mPAD_train_G(y,X2,Params,opts);

Z2Init = ones(size(X2))*0;
[X1,X2,Info] = mPAD_infer_X1X2_basis_func(y,Z2Init,Params,optsCG);
[X1FT,X2FT,ParamsFT] = mPAD_train_G_basis_func(y,Info.Z2,Params,opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Synthesis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TSamp = T; % number of samples of stimulus to generate

Dims = packDimsMPAD(D,K,TSamp); % Dimensionality of the sampled signal

[ySampICA,X1SampICA,X2SampICA,ASampICA] = sampleMPAD(ParamsICA,Dims);
[ySampSFA,X1SampSFA,X2SampSFA,ASampSFA] = sampleMPAD(ParamsSFA,Dims);
[ySampFT,X1SampFT,X2SampFT,ASampFT] = sampleMPAD(ParamsFT,Dims);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Checking the statistics of the sampled signals to see whether
% they (roughly) match the original

%% First plot the analysis
%figAnal=plotYX2(y,XSFA,'analysis',fs);

%% Plot the synthesis
%figAnal=plotYX2(ySamp,X2Samp,'synthesis',fs);

% Plot various statistics - you might need a fairly long sampled
% signal for these statistics to settle down. You can check this by
% looking how much they vary over repreated samples

y = y/sqrt(var(y));
ySampFT = ySampFT/sqrt(var(ySampFT));

statistics1 = get_statistics(y);
statistics2 = get_statistics(ySampFT);
statistics3 = get_statistics(ySampSFA);

plot_statistics(statistics1,statistics3);
plot_statistics(statistics1,statistics2);

%compareStatistics(y,ySampFT,C',X1SampFT,X2,X2SampFT,A',ASampFT,fs)

figure
hold on
plot(A(10,:));
plot(ASampFT(:,10),'-r')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SAVE
savePath = [saveDir,File(1:2),'/'];

baseName = [File(1:2),'_D',num2str(D),'_LrnLen_synth_basis_func_'];

name1 = [baseName,typ];
wavwrite(0.95*ySampFT/max(abs(ySampFT)),fs,[savePath,name1,'.wav']);
save([savePath,name1,'.mat'],'ParamsFT')

name2 = [baseName,'noFT_SFA'];
wavwrite(0.95*ySampSFA/max(abs(ySampSFA)),fs,[savePath,name2,'.wav']);
save([savePath,name2,'.mat'],'ParamsSFA')

name3 = [baseName,'noFT_ICA'];
wavwrite(0.95*ySampICA/max(abs(ySampICA)),fs,[savePath,name3,'.wav']);
save([savePath,name3,'.mat'],'ParamsICA')

name4 = [baseName,'orig'];
wavwrite(0.95*y/max(abs(y)),fs,[savePath,name4,'.wav']);

name6 = [baseName,'reversed'];
wavwrite(0.95*y(end:-1:1)/max(abs(y)),fs,[savePath,name6,'.wav']);

name5 = [baseName,'noise'];
wavwrite(0.95*ySampNoise/max(abs(ySampNoise)),fs,[savePath,name5,'.wav']);

['! mplayer ',name1,'.wav']