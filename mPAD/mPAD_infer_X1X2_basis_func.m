function [X1,X2,Info] = mPAD_infer_X1X2_basis_func(y,Z2,Params,varargin)

% function [X1,X2,Params,Info] = mPAD_train_G_basis_func(y,Z2,Params,varargin)
%
% Infers the modulators and carriers in mPAD. Requires initialised
% estimates of the envelopes and the parameters.
%
% INPUT
% y = Obsevations [T,1] 
% Z2 = initial estimates for transformed envelope coefficients [T,K] 
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

[T,K] = size(Z2);
D = length(Var1);


z2 = Z2(:);
  
L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

tol = 7;
X2 = Z2_to_X2(Z2,Len2,Var2,tol);

for l=1:L

  Gold = G;  
  X2old = X2;

  tic;
  [z2, ObjCur, itCur] = minimize(z2,'getObj_mPAD_noG_basis_func', ...
				 numIts(l),y,Params,tol);
  timCur = toc;


  Z2 = reshape(z2,[T,K]);      
  X2 = Z2_to_X2(Z2,Len2,Var2,tol);
    
  % plot if desired
  if verbose==1
    plot_mPAD_train_G(y,X2,X2old,Params)  
  end
  
  % Store objective and iteration information
  Obj = [Obj;ObjCur];
  it = [it;itCur];
  tim = [tim;timCur];
  dX2 = sqrt(sum((X2(:)-X2old(:)).^2)/sum(X2(:).^2));
  
  % Display some information to the user
  str1 = ['Progress ',num2str(l),'/',num2str(L),];
  str2 = ['Obj ',num2str(ObjCur(end))];
  if length(ObjCur)>1
    str3 = ['dObj ',num2str(diff(ObjCur(end-1:end)))];
  else
    str3 = ['dObj 0'];
  end
  str4 = ['time ',num2str(ceil(timCur/6)/10),'mins'];
  str5 = ['total time ',num2str(ceil(sum(tim)/6)/10),'mins'];
  str7 = ['dX2 ',num2str(round(dX2*1000)/10),'%%'];

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str7,str_space])

end

fprintf('\n')

% Save information about the optimisation
Info.Obj = Obj;
Info.it = it;
Info.Z2 = Z2;

% Estimate X1
A = log(1+exp(X2(1:T,:)*G'+ones(T,1)*Mu2'));
[ObjA,Xfin,Pfin] = kalman_mPAD_FB(Params,y,A');

[X1,covX1] = getFBLDSOutput(Xfin,Pfin);
X1 = X1';

