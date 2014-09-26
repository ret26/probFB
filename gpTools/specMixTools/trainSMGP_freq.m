function [varx,lenx,freqs,vary,info] = trainSMGP_freq(y,K,varargin)
%trainSMGP_freq trains a gp with a zero mean and a spectral mixture kernel
%  Inputs:
%    y: regularly sample data
%    K: number of components in the mixture
%  Outputs:
%    varx, lenx, freqs, vary: kernel hyperparameters and noise variance
%    info: some information about the optimisation, for debugging purposes
%  Optional inputs:
%    varargin{1}.numIts: number of iterations
%    varargin{1}.progress_chunk: report progress after 'chunk' iterations
%    varargin{2:5}: initialised params, varx, lenx, freqx and vary
%  See also getObjSMGP.m, initSMParams.m
% 
% Written by Richard Turner
% Modified by Thang Bui
% Last modified: 9/2014

if nargin>2 && isfield(varargin{1},'numIts')
    numIts = varargin{1}.numIts;
else
    numIts = 100;
end

if nargin>2 && isfield(varargin{1},'progress_chunk')
    progress_chunk = varargin{1}.progress_chunk;
else
    progress_chunk = 100;
end

specy = abs(fft(y)).^2;
T = length(y);
% initialisation if user doesn't provide starting values
if nargin<5
    % picking peaks in the freq domain, need to set win size and overlap
    params = initSMParams(y,K,100,0.5);
    varx = exp(params(1:3:3*K));
    lenx = exp(params(2:3:3*K));
    freqs = exp(params(3:3:3*K));
    vary = exp(params(end));
else
    % if user has provided initial values then use them
    varx = varargin{2};
    lenx = varargin{3};
    freqs = varargin{4};
    vary = varargin{5};
    params1 = [varx(:) lenx(:) freqs(:)]';
    params1 = params1(:);
    params = [params1; vary];
end
L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];



for l=1:L
    % save old parameters
    lenxOld = lenx;
    varxOld = varx;
    freqsOld = freqs;
    varyOld = vary;
    
    % TODO:
    
    % run conjugate gradient update
    tic;
    [params, ObjCur, itCur] = minimize(params,'getObjSMGP', ...
        numIts(l),specy);
    varx = exp(params(1:3:3*K));
    lenx = exp(params(2:3:3*K));
    freqs = exp(params(3:3:3*K));
    vary = exp(params(end));
    
    
    timCur = toc;
    % Store objective and iteration information
    Obj = [Obj;ObjCur];
    it = [it;itCur];
    tim = [tim;timCur];
    dvary = vary-varyOld;
    dlenx = lenx-lenxOld;
    dvarx = varx-varxOld;
    dfreqs = freqs - freqsOld;
    
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
    str6 = ['vary ',num2str(vary)];
    str7 = ['varx ',num2str(varx)];
    str8 = ['lenx ',num2str(lenx)];
    str9 = ['freqs ',num2str(freqs)];
    
    str_space = '   ';
    
    fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,...
        str_space,str5,str_space,str6,str_space,str7,str_space,str8,...
        str_space,str9]); 
end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
info.it = it;




end