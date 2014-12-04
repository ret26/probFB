function [f,df] = objFunctionSECosMix(theta,K,covfunc,X,Y,tau1,tau2,missingInd)

[fest,~] = predictSECosMix(theta,K,covfunc,X,Y,tau1,tau2,missingInd);
a = log10(abs(fft(fest)));
b = log10(abs(fft(Y)));
figure(10), 
plot((1:length(b)/2)/length(b),b(1:length(b)/2),'r')
hold on;
plot((1:length(a)/2)/length(a),a(1:length(a)/2),'b')
plot(exp(theta(2*K+1:end-1)),-2*ones(K,1),'or');
hold off;
drawnow
%exp(theta)

Ntrain = length(X);
noblks = Ntrain/tau1;
y1 = reshape(Y,tau1,noblks);
yKalman = num2cell(y1,1);


%% get model parameters and the derivatives
covhyp = reshape(theta(1:K*2),2,K)';
mu = exp(theta(K*2+1:end-1));
vary = exp(theta(end));
[Ct,Rt,At,Qt,dCt,dRt,dAt,dQt] = ...
    getParamsSECosMix(covfunc,covhyp,mu,vary,X,tau1,tau2,noblks,missingInd);

% initial condition
u1 = linspace(X(1),X(tau1),tau2)';
x0 = zeros(2*K*tau2,1);
P0 = zeros(2*K*tau2);
for k = 1:K
    Pk = feval(covfunc{:},covhyp(k,:),u1) + 1e-7*eye(tau2);
    P0(2*(k-1)*tau2+(1:tau2),2*(k-1)*tau2+(1:tau2)) = Pk;
    P0((2*k-1)*tau2+(1:tau2),(2*k-1)*tau2+(1:tau2)) = Pk;
end

dP0 = cell(length(theta),1);
for i = 1:length(theta)
    if i>K*2 % if freqs or observation noise
        dP0{i} = zeros(size(P0));
    else
        k = ceil(i/2);
        ind = i - (k-1)*2;
        dP0i = feval(covfunc{:},covhyp(k,:),u1,[],ind);
        dP0k = zeros(size(P0));
        dP0k(2*(k-1)*tau2+(1:tau2),2*(k-1)*tau2+(1:tau2)) = dP0i;
        dP0k((2*k-1)*tau2+(1:tau2),(2*k-1)*tau2+(1:tau2)) = dP0i;
        dP0{i} = dP0k;
    end
end

smoothing = 1;
verbose = 0;
[f,df] = kalmanVar(At,Ct,Qt,Rt,x0,P0,yKalman,verbose,...
    1-smoothing,dAt,dCt,dQt,dRt,dP0);
f = -f/(noblks*tau1);
df = -df/(noblks*tau1);
end

