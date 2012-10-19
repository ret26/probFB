function test_suite = test_mPAD_train_G
  initTestSuite;

% Tests: 
%
% function [X1,X2,Params] = mPAD_train_G(y,X2,Params,opts)
%

function test_runs

% quick test to make sure any modifications haven't broken the code

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings

seed = 1; 
randn('state',seed)
rand('state',seed)

T = 20;
FSamp = 8000;
D = 3;
K = 2;
Dims = packDimsMPAD(D,K,T);

% Carrier Processes
mVar1 = ones(1,D);
CF = [400,750,1000]/FSamp;
DF = [40,50,80]/FSamp;

% Modulator Processes
Len2 = [5,10]; % time constant in samples of modulators
Mu2 = randn(D,1);
Var2 = 9*ones(K,1);

G = [1,0;
     0,1;
     1,1];

vary = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DATA

[om,Lam1,Var1] = freq2probSpec(CF,DF,mVar1);

Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary,om);
[y,X1,X2,Amp] = sampleMPAD(Params,Dims);

opts.numIts = 6;
opts.alp = 1;
opts.progress_chunk = 3;
[X1,X2,Params,Info] = mPAD_train_G(y,X2,Params,opts);


function test_initialised_at_true_G_and_X2_stay_put

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings
dispFigs=1;

opts.verbose = dispFigs;
seed = 1; 
randn('state',seed)
rand('state',seed)

T = 1000;
FSamp = 8000;
D = 3;
K = 2;
Dims = packDimsMPAD(D,K,T);

% Carrier Processes
mVar1 = ones(1,D);
CF = [600,750,1000]/FSamp;
DF = [60,100,200]/FSamp;

% Modulator Processes
Len2 = [20,20]; % time constant in samples of modulators
Mu2 = -ones(D,1);
Var2 = 9*ones(K,1);

G = [1,0;
     0,1;
     1,1];

vary = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DATA

[om,Lam1,Var1] = freq2probSpec(CF,DF,mVar1);

Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary,om);
[y,X1,X2,Amp] = sampleMPAD(Params,Dims);

if dispFigs==1

  figure
  for d=1:D
    subplot(D,1,d)
    hold on
    Ychan = real(X1).*Amp;
    plot(Ychan(:,d),'-k')
    plot(Amp(:,d),'-r','linewidth',2)
  end
end

opts.numIts = 200; 
[X1Est,X2Est,ParamsEst,Info] = mPAD_train_G(y,X2,Params,opts);

GEst = ParamsEst.G;
AEst = log(1+exp(X2Est(1:T,:)*GEst'+ones(T,1)*Mu2'));

if dispFigs==1

  figure
  plot(Info.Obj)
  
  figure
  for d=1:D
    subplot(D,1,d)
    hold on
    YchanEst = real(X1Est).*AEst;
    
    plot(YchanEst(:,d),'-k')
    plot(Amp(:,d),'-b')
    plot(AEst(:,d),'-r','linewidth',2)  
  end
  
  figure
  subplot(1,2,1)
  imagesc(Params.G)

  subplot(1,2,2)
  imagesc(ParamsEst.G)

end

GError = sqrt(sum((G(:)-GEst(:)).^2))/sqrt(sum(G(:).^2));

if dispFigs==1
  disp(['Weights have moved by ',num2str(round(GError*100)),'%'])
end

assertTrue(GError<0.2) % weights move by less than 20%
%tol = 1e-1;
%assertVectorsAlmostEqual(G,GEst,'absolute',tol,0)


function test_initialised_at_random_G_and_true_X2

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward model Settings
dispFigs=1;
opts.verbose = dispFigs;

seed = 1; 
randn('state',seed)
rand('state',seed)

T = 1000;
FSamp = 8000;
D = 3;
K = 2;
Dims = packDimsMPAD(D,K,T);

% Carrier Processes
mVar1 = ones(1,D);
CF = [550,750,1000]/FSamp;
DF = [60,50,100]/FSamp;

% Modulator Processes
Len2 = [20,20]; % time constant in samples of modulators
Mu2 = -ones(D,1);
Var2 = 9*ones(K,1);

G = [1,0;
     0,1;
     1,1];

vary = 0.001;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GENERATE DATA

[om,Lam1,Var1] = freq2probSpec(CF,DF,mVar1);

Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary,om);
[y,X1,X2,Amp] = sampleMPAD(Params,Dims);

if dispFigs==1

  figure
  for d=1:D
    subplot(D,1,d)
    hold on
    Ychan = real(X1).*Amp;
    plot(Ychan(:,d),'-k')
    plot(Amp(:,d),'-r','linewidth',2)
  end
end

GScale = sqrt(sum(G.^2,1));
GInit = randn(D,K);
scale = diag(sqrt(diag(GInit'*GInit)));
GInit = (GInit/scale)*diag(GScale);
ParamsInit = Params;
ParamsInit.G = GInit;

opts.numIts = 200;
[X1Est,X2Est,ParamsEst,Info] = mPAD_train_G(y,X2,ParamsInit,opts);

GEst = ParamsEst.G;
AEst = log(1+exp(X2Est(1:T,:)*GEst'+ones(T,1)*Mu2'));

if dispFigs==1

  figure
  plot(Info.Obj)
  
  figure
  for d=1:D
    subplot(D,1,d)
    hold on
    YchanEst = real(X1Est).*AEst;
    
    plot(YchanEst(:,d),'-k')
    plot(Amp(:,d),'-b')
    plot(AEst(:,d),'-r','linewidth',2)  
  end
  
  figure
  subplot(1,2,1)
  imagesc(Params.G)

  subplot(1,2,2)
  imagesc(ParamsEst.G)

end

GError = sqrt(sum((G(:)-GEst(:)).^2))/sqrt(sum(G(:).^2));

disp(['Weights have moved by ',num2str(round(GError*100)),'%'])

assertTrue(GError<0.2) % weights move by less than 20%
%tol = 1e-1;
%assertVectorsAlmostEqual(G,GEst,'absolute',tol,0)


% function test_initialised_at_random_G_and_random_X2

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % Forward model Settings
% dispFigs=1;
% opts.verbose = dispFigs;
% seed = 1; 
% randn('state',seed)
% rand('state',seed)

% T = 4000;
% FSamp = 8000;
% D = 3;
% K = 3;
% Dims = packDimsMPAD(D,K,T);

% % Carrier Processes
% mVar1 = ones(1,D);
% CF = [1000,2000,3000]/FSamp;
% DF = [100,100,100]/FSamp;

% % Modulator Processes
% Len2 = [20,20,20]; % time constant in samples of modulators
% Mu2 = -ones(D,1);
% Var2 = 9*ones(K,1);

% G = [1,0,0;
%      0,1,1;
%      1,1,0];

% vary = 0.001;

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % GENERATE DATA

% [om,Lam1,Var1] = freq2probSpec(CF,DF,mVar1);

% Params = packParamsMPAD(G,Lam1,Var1,Len2,Var2,Mu2,vary,om);
% [y,X1,X2,Amp] = sampleMPAD(Params,Dims);

% if dispFigs==1

%   figure
%   for d=1:D
%     subplot(D,1,d)
%     hold on
%     Ychan = real(X1).*Amp;
%     plot(Ychan(:,d),'-k')
%     plot(Amp(:,d),'-r','linewidth',2)
%   end
% end

% GScale = sqrt(sum(G.^2,1));
% GInit = randn(D,K);
% scale = diag(sqrt(diag(GInit'*GInit)));
% GInit = (GInit/scale)*diag(GScale);
% ParamsInit = Params;
% ParamsInit.G = GInit;


% X2Init = 0*ones(T,K);

% [ObjA,Xfin,Pfin] = kalman_mPAD_FB(Params,y,ones(D,T));
% [X1,covX1] = getFBLDSOutput(Xfin,Pfin);
% AInit = abs(X1)';

% X2Init = log(exp(AInit)-1)-ones(T,1)*Params.Mu2';

% opts.numIts = 800;
% [X1Est,X2Est,ParamsEst,Info] = mPAD_train_G(y,X2Init,ParamsInit,opts);

% opts.numIts = 100;
% [X1Est2,X2Est2,ParamsEst2,Info2] = mPAD_train_G(y,X2,Params,opts);

% GEst = ParamsEst.G;
% AEst = log(1+exp(X2Est(1:T,:)*GEst'+ones(T,1)*Mu2'));

% if dispFigs==1

%   figure
%   plot(Info.Obj)
  
%   figure
%   for d=1:D
%     subplot(D,1,d)
%     hold on
%     YchanEst = real(X1Est).*AEst;
    
%     plot(YchanEst(:,d),'-k')
%     plot(Amp(:,d),'-b')
%     plot(AEst(:,d),'-r','linewidth',2)  
%   end
  
%   figure
%   subplot(1,2,1)
%   imagesc(Params.G)

%   subplot(1,2,2)
%   imagesc(ParamsEst.G)

% end

% GError = sqrt(sum((G(:)-GEst(:)).^2))/sqrt(sum(G(:).^2));

% disp(['Weights have moved by ',num2str(round(GError*100)),'%'])
% keyboard
% assertTrue(GError<0.2) % weights move by less than 20%
% %tol = 1e-1;
% %assertVectorsAlmostEqual(G,GEst,'absolute',tol,0)
