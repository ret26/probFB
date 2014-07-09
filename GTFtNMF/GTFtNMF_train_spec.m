function [W,H,lamv,varv,omv,mnV,covV,info] = GTFtNMF_train_spec(y,W,H,lamv,varv,omv,vary,lenx,mux,varx,varargin)

% function [W,H,mnV,covV,info] = GTFtNMF_train_spec(y,H,W,lamv,varv,omv,vary,lenx,mux,varx,varargin)
%
% Carries out learning by joint optimisation of the likelihood for the
% GTF-tNMF model. Unlike GTFtNMF this also learns the peripheral
% representation: that is varv, omv and lamv.
%
% y_t = sum_d Re(v_{d,t}) + sigy_t*\eta_t
% v_{d,t} = amp_{d,t}/amp_{d,t-1}*lamv_d*exp(i omv_d)*v_{d,t-1} 
%                              + amp_{d,t} \sigma_d \epsilon_{d,t}
% amp_{d,t} = sqrt(H*W/2)
% W_{k,d} = exp(logW_{k,d})/sum_d exp(logW_{k,d})
%
% Each of the temporal basis functions in NMF is given by: 
% H_{t,k} = exp(X_{t,k})  where X{1:T,k} ~ mux_k + GP(lenx_k,varx_k)
%
% The basis function representation of the GP used for the temporal
% basis functions is employed rather than the function-space view
% so that:
%
% bas_k \propto sqrt(varx_k)*exp(-1/(lenx_k^2) t^2)
% X_{1:T,k} = convolve(Z_{1:K,k},bas_k) 
% 
%
% INPUTS
% y = signal [T,1]
% W = initial spectral weight vector [K,D]
% H = initial temporal basis functions
% lamv = dynamical AR parameters [D,1]
% varv = dynamical noise parameters [D,1]
% omv = mean frequencies of the sinusoids [D,1]
% vary = observation noise [T,1]
% lenx = squared exponential length-scale [1,K]
% mux = steady state mean of the log temporal priors [1,K]
% varx = steady state variance of the log temporal priors [1,K]
%
% optional argument:
% Opts = structure of algorithmic options
%   tol = tolerance parameter which sets the extent of the Gaussian
%         basis functions
%   numIts = number of iterations for each M-Step [numEMIts,1]
%
% OUTPUTS
% W = learned spectral weight vector [K,D]
% H = learned temporal basis functions
% lamv = learned dynamical AR parameters [D,1]
% varv = learned dynamical noise parameters [D,1]
% omv = learned mean frequencies of the sinusoids [D,1]
% mnV = posterior mean over the 
%
% info = structure containing information about the optimisation including
%    Obj = objective function
%    it = number of completed iteration steps per batch
%    tim = time taken for each EM-iteration

if nargin>10 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 100;
end

if nargin>10 & isfield(varargin{1},'progress_chunk')
  progress_chunk = varargin{1}.progress_chunk;
else
  progress_chunk = 20;
end

if nargin>10 & isfield(varargin{1},'tol')
  tol = varargin{1}.tol;
else
  tol = 11;
end

disp('GTFtNMF: converting to basis function representation')

[K,D] = size(W);
[T,K] = size(H);
Z = zeros(T,K);

L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

for k=1:K
  % convert from 
  logH = log(H(:,k))-mux(k);
  [Z(:,k),infoTrans] = fitSEGP_BF(logH,lenx(k),varx(k));

  % logHfit = basis2SEGP(Z(:,k),lenx(k),varx(k),9);
  % figure
  % hold on
  % plot(logH,'-k')
  % plot(logHfit,'-b')
  % keyboard
end

% Initial amplitudes
A = (1/2*H*W).^(1/2);

% keeping tabs on progress
lik = [];

zlogWspec =  [Z(:);log(W(:));log(varv(:));log(lamv(:)./(1-lamv(:)));log((cos(omv(:))+1)./(1-cos(omv(:))))];

Obj = []; itTot = []; tim = []; lik = [];

for it=1:L
  
  % save old parameters
  Hold = H;
  Wold = W;
  Aold = A;
  varvold = varv;
  omvold = omv;
  lamvold = lamv;
  
  tic;
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % delta = 1e-6;
  % d=checkgrad('getObj_GTFtNMF_FB_hier_spec',zlogWspec,delta,y,vary, ...
  % 	    lenx,mux,varx,1,tol)
  % keyboard
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  % run conjugate gradient update
  tic;

  [zlogWspec, ObjCur, itCur] = minimize(zlogWspec,'getObj_GTFtNMF_FB_hier_spec', ...
				    numIts(it),y,vary,lenx,mux,varx,1,tol);
  timCur = toc;
%  keyboard
  
  % Pull out H, W and filter bank parameters from zlogW
  Z = reshape(zlogWspec(1:K*T),[T,K]);
  H = zeros([T,K]);
  
  for k=1:K
    tau = ceil(lenx(k)/sqrt(2)*tol);
    bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height to
						   % get correct variance
    bas = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
    
    H(:,k) =  exp(conv(Z(:,k),bas,'same')+mux(k));
  end
  
  W = reshape(exp(zlogWspec(T*K+1:T*K+K*D)),[K,D]);
  W = diag(1./sum(W,2))*W; %normalised weights

  q1 = zlogWspec(1+T*K+K*D:T*K+K*D+D);
  varv = exp(q1);
  
  q2 = zlogWspec(1+T*K+K*D+D:T*K+K*D+2*D);
  lamv = 1./(1+exp(-q2));

  q3 = zlogWspec(1+T*K+K*D+2*D:T*K+K*D+3*D);
  cosomv = (1-exp(-q3))./(1+exp(-q3));
  omv = angle(cosomv + i* sqrt(1-cosomv.^2));
  
  % Update amplitudes
  A = (1/2*H*W).^(1/2);

  itTot = [itTot;itCur];
  tim = [tim;timCur];
  Obj = [Obj;ObjCur];
  lik = [lik;-ObjCur(end)];
  
  dW = sqrt(sum((W(:)-Wold(:)).^2)/sum(W(:).^2));
  dH = sqrt(sum((H(:)-Hold(:)).^2)/sum(H(:).^2));
  dA = sqrt(sum((A(:)-Aold(:)).^2)/sum(A(:).^2));
  dlamv = sqrt(sum((lamv(:)-lamvold(:)).^2)/sum(lamv(:).^2));
  domv = sqrt(sum((omv(:)-omvold(:)).^2)/sum(omv(:).^2));
  dvarv = sqrt(sum((varv(:)-varvold(:)).^2)/sum(varv(:).^2));
  
  % Display some information to the user

  str1 = ['Progress ',num2str(it),'/',num2str(L)];
  str2 = ['lik ',num2str(-ObjCur(end))];

  if length(ObjCur)>1
    str3 = ['dObj ',num2str(diff(ObjCur(end-1:end)))];
  else
    str3 = ['dObj ','only one value'];
  end
  str4 = ['time ',num2str(ceil(timCur/6)/10),'mins'];
  str5 = ['total time ',num2str(ceil(sum(tim)/6)/10),'mins'];
  str6 = ['dH ',num2str(round(dH*1000)/10),'%%'];
  str7 = ['dW ',num2str(round(dW*1000)/10),'%%'];
  str8 = ['dA ',num2str(round(dA*1000)/10),'%%'];
  str9 = ['dlamv ',num2str(round(dlamv*1000)/10),'%%'];
  str10 = ['domv ',num2str(round(domv*1000)/10),'%%'];
  str11 = ['dvarv ',num2str(round(dvarv*1000)/10),'%%'];

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7,str_space,str8,str_space,str9,str_space,str10,str_space,str11])

end

fprintf('\n')

% Output posterior over Z
[likCur,Xfin,Pfin] = kalman_GTFtNMF_FB_hier(y,A',lamv,varv,omv,vary);
[mnV,covV] = getFBLDSOutput(Xfin,Pfin);

% Save information about the optimisation
info.lik = lik;
info.Obj = -Obj;
info.it = itTot;
info.tim = tim;
info.Z = Z;

