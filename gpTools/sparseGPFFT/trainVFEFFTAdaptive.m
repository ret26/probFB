function [params,info] = trainVFEFFTAdaptive(y,kernel,M,varargin)

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
    progress_chunk = 20;
end


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


% add some phantom datapoints to the start and end
T = length(y);
Textra = 2*round(exp(params(2)));
yextra = [zeros(Textra,1);y(:);zeros(Textra,1)];
T1 = T+2*Textra;

for l=1:L
    kernelspec = getGPSEisoSpec(exp(params(1)),exp(params(2)),T1);
    vary = exp(params(end));
    yextra = denoiseVFEFFT(yextra,kernelspec,vary,M);
    yextra(Textra+1:end-Textra) = y;
    figure(20), plot(1:T1,yextra,'r',Textra+1:T1-Textra,y,'g');
    specy = abs(fft(yextra)).^2;
    % run conjugate gradient update
    tic;
    [params, ObjCur, itCur] = minimize(params,'getObjVFEFFT', ...
        numIts(l),specy,kernel,M);
    
    timCur = toc;
    % Store objective and iteration information
    Obj = [Obj;ObjCur];
    it = [it;itCur];
    tim = [tim;timCur];
    %{
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
    %}
end

fprintf('\n')

% Save information about the optimisation
info.Obj = Obj;
info.it = it;

end