function [theta_end,nlml] = trainSE(theta_init,covfunc,...
    Xtrain,Ytrain,tau1,tau2,missingInd,noEvals)

Ntrain = length(Xtrain);
noblks = Ntrain/tau1;
y1 = reshape(Ytrain,tau1,noblks);
yKalman = num2cell(y1,1);
fname = 'objFunctionSE';

%test derivativevest(ind) = diag(C1*Pfint{i}*C1' + R1);
% d = checkgrad(fname,theta_init,1e-6,covfunc,...
%     Xtrain,yKalman,tau1,tau2,missingInd)
% keyboard

opt.length = -noEvals;
opt.method = 'BFGS';
opt.verbosity = 2;
[theta_end,lik,i] = minimize_new(theta_init,fname,opt,covfunc,...
    Xtrain,yKalman,tau1,tau2,missingInd);
nlml = lik(end);
end