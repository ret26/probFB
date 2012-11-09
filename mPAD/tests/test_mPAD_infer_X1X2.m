function test_suite = test_mPAD_infer_X1X2
  initTestSuite;

% Tests: 
%
% [X1,X2,Info] = mPAD_infer_X1X2(y,X2,Params,varargin)
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
opts.progress_chunk = 3;
[X1,X2,Info] = mPAD_infer_X1X2(y,X2,Params,opts);


function test_initialised_at_true_X2_stay_put

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
[X1Est,X2Est,Info] = mPAD_infer_X1X2(y,X2,Params,opts);

AEst = log(1+exp(X2Est(1:T,:)*G'+ones(T,1)*Mu2'));

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
  
end

X2Error = sqrt(sum(sum((X2(1:T,:)-X2Est(1:T,:)).^2)))/sqrt(sum(sum(X2(1:T,:).^2)));

if dispFigs==1
  disp(['Envelopes have moved by ',num2str(round(X2Error*100)),'%'])
end

assertTrue(X2Error<0.2) % envelopes move by less than 20%
%tol = 1e-1;
%assertVectorsAlmostEqual(G,GEst,'absolute',tol,0)


function test_initialised_at_random_X2_and_move_toward_true

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

X2Init = X2+randn(size(X2))/1;

opts.numIts = 200; 
[X1Est,X2Est,Info] = mPAD_infer_X1X2(y,X2Init,Params,opts);

AEst = log(1+exp(X2Est(1:T,:)*G'+ones(T,1)*Mu2'));

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
  
end

X2Error1 = sqrt(sum(sum((X2(1:T,:)-X2Est(1:T,:)).^2)))/sqrt(sum(sum(X2Est(1:T,:).^2)));

X2Error2 = sqrt(sum(sum((X2Init(1:T,:)-X2Est(1:T,:)).^2)))/sqrt(sum(sum(X2Est(1:T,:).^2)));

X2Error0 = sqrt(sum(sum((X2Init(1:T,:)-X2(1:T,:)).^2)))/sqrt(sum(sum(X2Est(1:T,:).^2)));


if dispFigs==1
  disp(['Envelopes have moved by ',num2str(round(X2Error2*100)),'%'])
  disp(['Error on initialisation ',num2str(round(X2Error0*100)),'%'])
  disp(['Error after training ',num2str(round(X2Error1*100)),'%'])
end

assertTrue(X2Error1<X2Error0) % envelopes move closer to the true
                              % values - not a very strong test 


