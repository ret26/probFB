function [W,H,mnV,covV,info] = GTFtNMF_EM(y,W,H,lamv,varv,omv,vary,lenx,mux,varx,varargin)

% function [W,H,mnV,covV,info] = GTFtNMF_EM(y,H,W,lamv,varv,omv,vary,lenx,mux,varx,varargin)
%
% Carries out expectation maximisation learning of the GTF-tNMF
% model
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
%   numEMIts = number of iterations of EM
%   numIts = number of iterations for each M-Step [numEMIts,1]
%
% OUTPUTS
% W = initial spectral weight vector [K,D]
% H = initial temporal basis functions
% mnV = posterior mean over the 
% covV = posterior covariance of the time-frequency coefficients
% info = structure containing information about the optimisation including
%    Obj = objective function
%    it = number of completed iteration steps per batch
%    tim = time taken for each EM-iteration

if nargin>10 & isfield(varargin{1},'numEMIts')
  numEMIts = varargin{1}.numEMIts;
else
  numEMIts = 50;
end

if nargin>10 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = ones(numEMIts,1)*50;
end

if nargin>10 & isfield(varargin{1},'tol')
  tol = varargin{1}.tol;
else
  tol = 11;
end

disp('GTFtNMF_EM: converting to basis function representation')

[K,D] = size(W);
[T,K] = size(H);
Z = zeros(T,K);

for k=1:K
  % convert from 
  logH = log(H(:,k))-mux(k);
  [Z(:,k),infoTrans] = fitSEGP_BF(logH,lenx(k),varx(k));
  H(:,k) = exp(basis2SEGP(Z(:,k),lenx(k),varx(k),tol)+mux(k));
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

zlogW = [Z(:);log(W(:))];
Obj = []; itTot = []; tim = [];

for it=1:numEMIts
  
  % save old parameters
  Hold = H;
  Wold = W;
  Aold = A;
  tic;

  % E-Step
  [likCur,VV,RVV] = kalman_GTFtNMF_FB_hier_EM(y,A',lamv,varv,omv,vary);
  likCur = (likCur-1/2*sum(Z(:).^2))/T; % add the prior on and
                                        % scale by # data points

% after the first iteration the likelihood above is equal to
%  likCur2 = -getObj_GTFtNMF_FB_hier(zlogW,y,lamv,varv,omv, ...
%				    vary,lenx,mux,varx,tol);

  % M-Step
  [zlogW, ObjCur, itCur] = minimize(zlogW,'getObj_GTFtNMF_FB_hier_EM', ...
				    numIts(it),VV',RVV',lamv,varv, ...
				    lenx,mux,varx,1,tol);
  timCur = toc;
  
  
  
  % Pull out H and W from zlogW
  Z = reshape(zlogW(1:K*T),[T,K]);
  H = zeros([T,K]);

  % store information about progress -- remember to add prior to likelihood
  lik = [lik;likCur];

  for k=1:K
    tau = ceil(lenx(k)/sqrt(2)*tol);
    bh = sqrt(varx(k)*sqrt(2)/(lenx(k)*sqrt(pi))); % basis function height to
						   % get correct variance
    bas = bh*exp(-1/(lenx(k)^2)*([-tau:1:tau]').^2);
    
    H(:,k) =  exp(conv(Z(:,k),bas,'same')+mux(k));
  end
  
  W = reshape(exp(zlogW(K*T+1:end)),[K,D]);
  W = diag(1./sum(W,2))*W;
  
  % Update amplitudes
  A = (1/2*H*W).^(1/2);

  itTot = [itTot;itCur];
  tim = [tim;timCur];
  Obj = [Obj;ObjCur];
  
  dW = sqrt(sum((W(:)-Wold(:)).^2)/sum(W(:).^2));
  dH = sqrt(sum((H(:)-Hold(:)).^2)/sum(H(:).^2));
  dA = sqrt(sum((A(:)-Aold(:)).^2)/sum(A(:).^2));
  
  
  % Display some information to the user
%keyboard
  str1 = ['Progress ',num2str(it),'/',num2str(numEMIts),];
  str2 = ['lik ',num2str(likCur)];

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

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7,str_space,str8])

end

fprintf('\n')

% Output posterior over Z
[likCur,Xfin,Pfin] = kalman_GTFtNMF_FB_hier(y,A',lamv,varv,omv,vary);
[mnV,covV] = getFBLDSOutput(Xfin,Pfin);
lik = [lik;(likCur-1/2*sum(Z(:).^2))/T];

% Save information about the optimisation
info.Obj = Obj;
info.it = itTot;
info.tim = tim;
info.lik = lik;
info.Z = Z;
