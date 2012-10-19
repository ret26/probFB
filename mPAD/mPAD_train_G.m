function [X1,X2,Params,Info] = mPAD_train_G(y,X2,Params,varargin)

% function [X1,X2,Params,Info] = mPAD_train_G(y,X2,Params,varargin)
%
% Trains the weight parameters in mPAD - only alters the directions
% and not the magnitudes. Requires initialised estimates of the
% envelopes and parameters.
%
% INPUT
% y = Obsevations [T,1] 
% X2 = Transformed envelopes [Tx,K] 
% Params = structure containing the following parameters
%   G = modulator weights [D,K]
%   Lam1 = Carrier dynamics [D,2]
%   Var1 = Carrier noise [1,D]
%   Len2 = Modulator Length Scales [K,1] (not required)
%   Var2 = Modulator Noise [1,K]
%   vary = observation noise  
% opts = structure of optional options
%   numIts = number of iterations
%   progress_chuck = number of iterations for each conjugate
%                    gradient run
%   verbose = if verbose==1 plots graphs each iteration
%
% OUTPUT
% X1 = carriers (complex) [T,1]
% X2 = Transformed envelopes [Tx,K]
% Params = structure of parameters, every field equal to the input
%          except for G which will have been adapted to match the
%          signal
% Info = structure containing information about the optimisation
%        Obj = objectives
%        it = number of completed iterations

if nargin>3 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 200;
end

if nargin>3 & isfield(varargin{1},'progress_chunk')
  progress_chunk = varargin{1}.progress_chunk;
else
  progress_chunk = 10;
end

if nargin>3 & isfield(varargin{1},'verbose')
  verbose = varargin{1}.verbose;
else
  verbose = 0;
end

[G,Lam1,Var1,Len2,Var2,Mu2,vary,om] = ...
    unpackParamsMPAD(Params);

T = length(y);
[Tx,K] = size(X2);
D = length(Var1);

% If the user just intialises the first T time-steps of the
% envelopes, then we need to augment the X2 vector

if T==Tx
  tau = 5*max(Len2); % offset added to avoid wrap-around effects 
  Tx = 2^ceil(log2(T+tau)); % make the duration a power of 2 so fft is faster
  X2 = addRamp(X2,Len2); % ramp the onset and offset to smooth the
                         % circular boundary conditions
  X2 = [X2;zeros(Tx-T,K)]; % add the extra-zero padding
end

% compute the spectra of the GPs
fftCov = getGPSESpec(Len2,Tx);

% store the scales of the weights
GScale = sqrt(sum(G.^2,1));

% minimise the objective function
X2G = [X2(:);G(:)];

L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];


for l=1:L

  Gold = G; 
  X2old = X2;

  tic;
  [X2G, ObjCur, itCur] = minimize(X2G,'getObj_mPAD_fixG', ...
				  numIts(l),y,fftCov, Params,GScale);
  timCur = toc;


  % Pull out G and X2 from X2G
  G = reshape(X2G(Tx*K+1:Tx*K+K*D),[D,K]);
  scale = diag(sqrt(diag(G'*G)));
  G = (G/scale)*diag(GScale);
  Params.G = G;

  % need to figure out whether to return X2(1:T,:) or X2(1:Tx,:)??
  
  X2 = reshape(X2G(1:Tx*K,1),[Tx,K]);    
  % dX2new = X2-X2old;
  % X2new = -dX2new/5+X2;
  % X2G = [X2new(:);G(:)];
    
  % plot if desired
  if verbose==1
    plot_mPAD_train_G(y,X2,X2old,Params)  
  end
  
  % Store objective and iteration information
  Obj = [Obj;ObjCur];
  it = [it;itCur];
  tim = [tim;timCur];
  dG = sqrt(sum((G(:)-Gold(:)).^2)/sum(G(:).^2));
  dX2 = sqrt(sum((X2(:)-X2old(:)).^2)/sum(X2(:).^2));
  
  % Display some information to the user
  str1 = ['Progress ',num2str(l),'/',num2str(L),];
  str2 = ['Obj ',num2str(ObjCur(end))];
  str3 = ['dObj ',num2str(diff(ObjCur(end-1:end)))];
  str4 = ['time ',num2str(ceil(timCur/6)/10),'mins'];
  str5 = ['total time ',num2str(ceil(sum(tim)/6)/10),'mins'];
  str6 = ['dG ',num2str(round(dG*1000)/10),'%%'];
  str7 = ['dX2 ',num2str(round(dX2*1000)/10),'%%'];

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str7,str_space])

end

fprintf('\n')

% Save information about the optimisation
Info.Obj = Obj;
Info.it = it;


% Estimate X1
A = log(1+exp(X2(1:T,:)*G'+ones(T,1)*Mu2'));
[ObjA,Xfin,Pfin] = kalman_mPAD_FB(Params,y,A');

[X1,covX1] = getFBLDSOutput(Xfin,Pfin);
X1 = X1';