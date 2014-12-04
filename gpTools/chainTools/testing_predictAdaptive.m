function [fest,vest] = testing_predictAdaptive(theta,K,covfunc,Xtrain,Ytest,tau1,tau2,missingInd)

Ntrain = length(Xtrain);
noblks = Ntrain/tau1;
y1 = reshape(Ytest,tau1,noblks);
yKalman = num2cell(y1,1);

% get model parameters
covhyp = reshape(theta(1:K*2),2,K)';
mu = exp(theta(K*2+1:end-1));
vary = exp(theta(end));
[Ct,Rt,At,Qt] = ...
    getParamsSECosMix(covfunc,covhyp,mu,vary,Xtrain,tau1,tau2,noblks,missingInd);
% initial condition
u1 = linspace(Xtrain(1),Xtrain(tau1),tau2)';
x0 = zeros(2*K*tau2,1);
P0 = zeros(2*K*tau2);
for k = 1:K
    Pk = feval(covfunc{:},covhyp(k,:),u1) + 1e-7*eye(tau2);
    P0(2*(k-1)*tau2+(1:tau2),2*(k-1)*tau2+(1:tau2)) = Pk;
    P0((2*k-1)*tau2+(1:tau2),(2*k-1)*tau2+(1:tau2)) = Pk;
end
% keyboard
[Xfint,Pfint] = testing_kalmanAdaptive(At,Ct,Qt,Rt,x0,P0,yKalman,1,0);

u = linspace(1,tau1,tau2);
s = 1:tau1;
tin = 1e-7;
Ck = cell(K,1);
Rk = cell(K,1);
for k = 1:K
    Kss = feval(covfunc{:},covhyp(k,:),s') + tin*eye(tau1);
    Kuu = feval(covfunc{:},covhyp(k,:),u') + tin*eye(tau2);
    Ksu = feval(covfunc{:},covhyp(k,:),s',u');
    L0 = chol(Kuu);
    Ksu_invL0 = Ksu/L0;
    Ck{k} = Ksu_invL0/L0';
    Rk{k} = Kss - Ksu_invL0*Ksu_invL0';
end

t = noblks;
ind = tau1*(t-1) + (1:tau1);
Ct = [];
R1 = zeros(tau1);
for k = 1:K
    mukInd = mu(k)*ind(:);
    cosVec = cos(2*pi*mukInd);
    sinVec = sin(2*pi*mukInd);
    cosMat11 = cosVec(:,ones(1,tau1));
    sinMat11 = sinVec(:,ones(1,tau1));
    cosMat12 = cosVec(:,ones(1,tau2));
    sinMat12 = sinVec(:,ones(1,tau2));
    Ct = [Ct cosMat12.*Ck{k} sinMat12.*Ck{k}];
    R1 = R1 + cosMat11.*Rk{k}.*cosMat11' + sinMat11.*Rk{k}.*sinMat11'; 
end

fest = Ct*Xfint{t};
vest = diag(Ct*Pfint{t}*Ct' + R1);
end

