function [f,df] = vfeObjFunc(theta,covfunc,x,y,varargin)

%
[n,D] = size(x);
noParams1  = eval(feval(covfunc{:}));
optU = 0;
if length(theta) ~= noParams1 + 1
    optU = 1;
end

covParams = theta(1:noParams1);
sigma  = exp(theta(noParams1+1));
if optU
    noUParams = length(theta) - noParams1 - 1;
    m = noUParams/D;
    z = theta(noParams1+2:end);
    z = reshape(z,m,D);
else
    z = varargin{1};
    m = size(z,1);
end

tiny = 1e-7;

% compute the covariance matrix
Knm = feval(covfunc{:},covParams,x,z);
Kmn = Knm';
Kmm = feval(covfunc{:},covParams,z);
Knndiag = feval(covfunc{:},covParams,x,'diag');


% compute the objective function
sigma2 = sigma^2;
cterm = 1/2/sigma2;
f0 = -n/2*log(2*pi) - (n-m)*log(sigma) - cterm*(y'*y);
cholKmm = chol(Kmm+tiny*eye(m));
invUmm =  cholKmm\eye(m);
f1 = sum(log(diag(cholKmm)));
A1 = Knm'*Knm;
A = sigma2*Kmm + A1;
A = A + tiny*eye(size(A));
cholA = chol(A);
invUA =  cholA\eye(m);
f2 = -sum(log(diag(cholA)));
KmnY = Kmn*y;
Xb = solve_chol(cholA,KmnY);
f3 = cterm*KmnY'*Xb;
f4 = -cterm*sum(Knndiag);
KnmInvUmm = Knm*invUmm;
C = KnmInvUmm'*KnmInvUmm;
f5 = cterm*trace(C);
%f1+f2, f3, f4+f5
f = f0 + f1 + f2 + f3 + f4 + f5;

invKmm = invUmm*invUmm';
invA   = invUA*invUA';
Tmm = invKmm - sigma2*invA - Xb*Xb';
T1 = Tmm - invKmm*A1*invKmm/sigma2;
T2 = Tmm*Kmn + Xb*y';
% compute the derivatives
df = [];
for j = 1:length(theta)
    if j <= length(covParams) % kernel hyper
        dKmmj = feval(covfunc{:},covParams,z,[],j);
        dKnmj = feval(covfunc{:},covParams,x,z,j);
        dKnnj = feval(covfunc{:},covParams,x,'diag',j);
        dj = 0.5*trace(dKmmj*T1) - cterm*sum(dKnnj) + trace(T2*dKnmj)/sigma2;
        df = [df; dj];
    elseif j == length(covParams)+1 % noise variance
        dj1 = Xb'*Kmm*Xb;
        dj = -(n-m) + y'*y/sigma2 - 2*(f4+f5+f3) - dj1 - sigma2*trace(invA*Kmm);
        df = [df; dj];
    end
end

if optU
    % only support covSEisoU, covSEiso and covSEard atm
    if iscell(covfunc), covstr = covfunc{1}; else covstr = covfunc; end
    if ~ischar(covstr), covstr = func2str(covstr); end
    if ~strcmp(covstr,'covSEisoU') && ~strcmp(covstr,'covSEiso') ...
            && ~strcmp(covstr,'covSEard')
        error('dL/dz only supported for covSEisoU, covSEiso and covSEard!');
    end
    D1 = zeros(D,m);
    D2 = zeros(D,m);
    for d=1:D
        DKnm = ((x(:,d)*ones(1,m) - ones(n,1)*z(:,d)')/exp(2*theta(d))).*Knm;
        DKmm = -((ones(m,1)*z(:,d)' - z(:,d)*ones(1,m))/exp(2*theta(d))).*Kmm;
        D1(d,:) = sum(DKmm.*T1, 1);
        D2(d,:) = sum(DKnm.*T2', 1);
    end
    Du = D1 + D2/sigma2;
    Du = reshape(Du',m*D,1);
    df = [df; Du]; % TODO: check!
end

% 
f = -f;
df = -df;
% keyboard

end
