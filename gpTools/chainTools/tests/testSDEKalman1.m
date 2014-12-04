
close all; clear all; clc
%% Generate data
% Discretize x (measurement points)
x = (-50:0.2:50)';

% Test points (evaluation points)
xt = linspace(-60,60,1000)';

sigma2 = .1;

covfunc = {@covSEiso};
hyp.cov = log([5 1]);
Kxx = feval(covfunc{:},hyp.cov,x);
N = length(x);
y1 = chol(Kxx+1e-6*eye(N))'*randn(N,1);
y = y1 + sqrt(sigma2)*randn(N,1);


%% full GP using gpml
hyp1.cov = log([5 1]);
hyp1.lik = log(0.1);
hyp2 = minimize(hyp1,@gp,-100,@infExact,[],covfunc,@likGauss,x,y);

[mf,sf] = gp(hyp2,@infExact,[],covfunc,@likGauss,x,y,xt);
func1 = @() gp(hyp2,@infExact,[],covfunc,@likGauss,x,y,xt);
timeGPML = timeit(func1)
figure(1), plot(x,y,'+b',xt,mf,'or')


%% kalman GP using gpstuff 

% The likelihood model
lik = lik_gaussian('sigma2', 2);


% The GP covariance function
gpcf = gpcf_sexp('lengthScale', 5, ...
    'magnSigma2', 1, 'kalman_deg', 1);

% Define Gaussian process model using type 'KALMAN'
gp1 = gp_set('lik', lik, 'cf', gpcf, 'type', 'KALMAN');
% Hyperparameter optimization (state space model)
gp1 = gp_optim(gp1, x, y);

% Predict values at test inputs xt
[m1,v1] = gp_pred(gp1, x, y, 'xt', xt);
func2 = @() gp_pred(gp1, x, y, 'xt', xt);
timeKalman = timeit(func2)
%% full gp using gpstuff
% Hyperparameter optimization (full model)
gp2 = gp_set(gp1,'type','FULL');
gp_full = gp_optim(gp2, x, y);
% Predict values at test inputs xt
[m2,v2] = gp_pred(gp_full, x, y, 'xt', xt);
func3 = @() gp_pred(gp_full, x, y, 'xt', xt);
timeGPStuff = timeit(func3)

%% Compare against full GP solution (table)

% Show table with comparison to full GP results
fprintf('\n%12s | %8s | %8s | %8s \n', ...
    'Parameter','FULL GPML', 'FULL GPSTUFF', 'KALMAN')

fprintf('-----------------------------------\n')

fprintf('%12s | %8.4f | %8.4f | %8.4f \n', ...
    'magnSigma2',exp(hyp2.cov(2)),gp_full.cf{1}.magnSigma2,gp1.cf{1}.magnSigma2)
fprintf('%12s | %8.4f | %8.4f | %8.4f \n', ...
    'lengthScale',exp(hyp2.cov(1)),gp_full.cf{1}.lengthScale,gp1.cf{1}.lengthScale)
fprintf('%12s | %8.4f | %8.4f | %8.4f \n', ...
    'sigma2',exp(hyp2.lik),gp_full.lik.sigma2,gp1.lik.sigma2)

d1 = sum(abs(mf-m1))
d2 = sum(abs(mf-m2))

%% Show result

figure(2), 
plot(x,y,'+b',xt,m1,'ok')
figure(3),
plot(x,y,'+b',xt,m2,'om')

%legend('Input','GPML','GPSTUFF','KALMAN');

%xlabel('Input, x');
%ylabel('Output, y');


% %% Compare the full posteriors
% 
% % Return the posterior mean EFT and covariance COVFT of latent
% % variables:
% %     Eft =  E[f | xt,x,y,th]  = K_fy*(Kyy+s^2I)^(-1)*y
% %   Covft = Var[f | xt,x,y,th] = K_fy - K_fy*(Kyy+s^2I)^(-1)*K_yf.
% 
% figure(2); clf
% 
% % The full GP
% subplot(121)
% 
% % Calculate
% [Eft1,Covft] = gp_jpred(gp_full,x,y,xt);
% 
% % Plot
% imagesc(Covft)
% axis equal tight
% 
% % Title
% title('The full GP posterior covariance')
% 
% % The state space model ('KALMAN')
% subplot(122)
% 
% % Calculate
% [Eft1,Covft] = gp_jpred(gp,x,y,xt);
% 
% % Plot
% imagesc(Covft)
% axis equal tight
% 
% % Title
% title('The state space GP posterior covariance')
