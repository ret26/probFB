function [params,info] = trainSparseGP(y,kernel,M,varargin)

% kernel.name = SE or SM
% for SM, kernel.K specifies the number of components in the mixture


if nargin>3 && isfield(varargin{1},'numIts')
    numIts = varargin{1}.numIts;
else
    numIts = 100;
end

if nargin>3 && isfield(varargin{1},'progress_chunk')
    progress_chunk = varargin{1}.progress_chunk;
else
    progress_chunk = 100;
end

specy = abs(fft(y)).^2;
T = length(y);
% initialisation if user doesn't provide starting values
if nargin<5
    % TODO: init params depending on the kernel provided
    error('TODO: init params')
else
    % if user has provided initial values then use them
    params = varargin{2};
end
L = ceil(numIts/progress_chunk);
numIts = ones(L,1)*progress_chunk;

Obj = []; it = []; tim = [];

for l=1:L
    
    % run conjugate gradient update
    tic;
    [params, ObjCur, itCur] = minimize(params,'getObjSparseGP', ...
        numIts(l),specy,kernel,M);
    
    timCur = toc;
    % Store objective and iteration information
    Obj = [Obj;ObjCur];
    it = [it;itCur];
    tim = [tim;timCur];
    
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
    
    str_space = '   ';
    
    fprintf(['\n',str1,str_space,str2,str_space,str3,str_space,str4,...
        str_space,str5]); 
end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
info.it = it;

end