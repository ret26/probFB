function [y,aRec,info] = recon_FB_mag(y,aTar,spec,varargin)

% function [y,info] = recon_FB_mag(y,aTar,spec,varargin)
%
% reconstructs a signal from amplitudes taken from the subbands of
% the signal (equivalently spectrogram coefficients)
%
% INPUTS
% y = initial guess at the signal [T,1]
% aTar = target envelopes [T,D]
% spec = subband filter ffts [T,D]
%
% Optional algorithmic options:
% numIts = number of iterations per optimisation batch
% progress_chunk = number of optimisation batches 
% error_measure = 'a' (squared error in evelopes) 
%                 'loga' (log error in envelopes)
% 
% [y,aRec,info] = recon_FB_mag(y,aTar,spec,opts,lam,vary)
% This will impose a prior (regularisation) onto y
% 
% OUTPUTS
% y = reconstructed signal [T,1]
% aRec = reconstructed signal envelopes [T,D]
%   info = structure containing information about the optimisation including
%     Obj = objective function
%     it = number of completed iteration steps per batch
%
if nargin>3 & isfield(varargin{1},'numIts')
  numIts = varargin{1}.numIts;
else
  numIts = 1000;
end

if nargin>3 & isfield(varargin{1},'progress_chunk')
  progress_chunk = varargin{1}.progress_chunk;
else
  progress_chunk = 50;
end

if nargin>3 & isfield(varargin{1},'error_measure')
  flag = varargin{1}.error_measure;
else
  flag = 'a';
end

L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

for l=1:L
  
  % save old signal
  yold = y;
  
  % run conjugate gradient update
  tic;
  if nargin<5
    [y, ObjCur, itCur] = minimize(y,'get_obj_recon_FB_mag', ...
				      numIts(l),aTar,spec,flag);
  else
    lam = varargin{2};
    vary = varargin{3};
    [y, ObjCur, itCur] = minimize(y,'get_obj_recon_FB_mag', ...
				  numIts(l),aTar,spec,flag,lam,vary);
  end
    timCur = toc;

  % Store objective and iteration information
  Obj = [Obj;ObjCur];
  it = [it;itCur];
  tim = [tim;timCur];
  dy = sqrt(sum((y(:)-yold(:)).^2)/sum(y.^2));
  
  % Display some information to the user
  str1 = ['Progress ',num2str(l),'/',num2str(L),];
  str2 = ['Obj ',num2str(ObjCur(end))];

  if length(ObjCur)>1
    str3 = ['dObj ',num2str(diff(ObjCur(end-1:end)))];
  else
    str3 = ['dObj ','only one value'];
  end
  str4 = ['time ',num2str(ceil(timCur/6)/10),'mins'];
  str5 = ['total time ',num2str(ceil(sum(tim)/6)/10),'mins'];
  str6 = ['dy ',num2str(round(dy*1000)/10),'%%'];

  str_space = '   ';
  
  fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,str_space,str5,str_space,str6,str_space,str_space])

end

fprintf('\n')

% compute reconstructed signal's envelopes
[T,D] = size(aTar);
yFFT = fft(y);

aRec = zeros(T,D);

for d=1:D
  zCur = ifft(spec(:,d).*yFFT);
  aRec(:,d) = abs(zCur);
end

% Save information about the optimisation
info.Obj = Obj;
info.it = it;
