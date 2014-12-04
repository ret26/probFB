function [theta_end,nlml] = testing_trainAdaptive(theta_init,K,covfunc,...
    Xtrain,Ytrain,tau1,tau2,missingInd,noEvals,weights)
fname = 'testing_objFunctionAdaptive';

% %test derivative
% d = checkgrad(fname,theta_init,1e-6,K,covfunc,...
%     Xtrain,Ytrain,tau1,tau2,missingInd)
% keyboard

opt.length = -noEvals;
opt.method = 'CG';
opt.verbosity = 2;
[theta_end,lik,i] = minimize_new(theta_init,fname,opt,K,covfunc,...
    Xtrain,Ytrain,tau1,tau2,missingInd,weights);


nlml = lik(end);
end