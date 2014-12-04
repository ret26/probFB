function [fest,vest] = predictSELocal(theta,covfunc,Xtrain,Ytest,tau1,tau2,missingInd)

Ntrain = length(Xtrain);
noblks = Ntrain/tau1;
y1 = reshape(Ytest,tau1,noblks);
yKalman = num2cell(y1,1);

% get model parameters
vary = exp(theta(end));
[Ct,Rt,Qt] = ...
    getParamsSELocal(covfunc,theta(1:end-1),vary,Xtrain,tau1,tau2,noblks,missingInd);

% initial condition
u1 = linspace(Xtrain(1),Xtrain(tau1),tau2)';
x0 = zeros(tau2,1);
D = size(Xtrain,2);
alpha = eval(feval(covfunc{:}));
P0 = feval(covfunc{:},theta(1:alpha),u1);
% keyboard
[Xfint,Pfint] = kalmanVarLocal(Ct,Qt,Rt,x0,P0,yKalman,1,0);

fest = zeros(size(Ytest));
vest = zeros(size(Ytest));
u = linspace(1,tau1,tau2);
s = 1:tau1;
tin = 1e-6;
hyp.cov = theta(1:end-1);
Kss = feval(covfunc{:},hyp.cov,s') + tin*eye(length(s));
Kuu = feval(covfunc{:},hyp.cov,u') + tin*eye(length(u));
Ksu = feval(covfunc{:},hyp.cov,s',u');
C = Ksu/Kuu;
R = Kss - Ksu/Kuu*Ksu';
for i = 1:noblks
    ind = tau1*(i-1) + (1:tau1);
    fest(ind) = C*Xfint{i};
    vest(ind) = diag(C*Pfint{i}*C' + R);
end
end

