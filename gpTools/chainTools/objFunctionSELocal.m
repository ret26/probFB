function [f,df] = objFunctionSELocal(theta,covfunc,X,Y,tau1,tau2,missingInd)


%% get model parameters and the derivatives
noblks = length(Y);
vary = exp(theta(end));
[Ct,Rt,Qt,dCt,dRt,dQt] = ...
    getParamsSELocal(covfunc,theta(1:end-1),vary,X,tau1,tau2,noblks,missingInd);

% initial condition
u1 = linspace(X(1),X(tau1),tau2)';
x0 = zeros(tau2,1);
D = size(X,2);
alpha = eval(feval(covfunc{:}));
P0 = feval(covfunc{:},theta(1:alpha),u1) + 1e-7*eye(tau2);
dP0 = cell(alpha+1,1);
for i = 1:alpha
    dP0i = feval(covfunc{:},theta(1:alpha),u1,[],i);
    dP0{i} = dP0i;
end
dP0{end} = zeros(size(P0));
smoothing = 1;
verbose = 0;
[f,df] = kalmanVarLocal(Ct,Qt,Rt,x0,P0,Y,verbose,...
    1-smoothing,dCt,dQt,dRt,dP0);
f = -f/(noblks*tau1);
df = -df/(noblks*tau1);
end

