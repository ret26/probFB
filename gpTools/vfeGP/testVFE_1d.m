% test vfe code
clear all, close all

%% generate data using gpml package
rng(234);
disp('generating data ...');
m = 10;
covfunc = {@covSEiso}; ell = 0.3; sf = 1; hyp.cov = log([ell; sf]);
likfunc = @likGauss; sn = 0.1; hyp.lik = log(sn);

n = 100;
x = 5*randn(n, 1);
K = feval(covfunc{:},hyp.cov,x) + 1e-6*eye(n);
y = chol(K)'*randn(n, 1) + exp(hyp.lik)*randn(n, 1);


%% FULL GP
disp('running full gp ...');
z = linspace(-15,15,1000)';
%z = sort(x);
hyp1 = minimize(hyp, @gp, -1000, @infExact, [], covfunc, likfunc, x, y);
[m1,s21] = gp(hyp1, @infExact, [], covfunc, likfunc, x, y, z);
f1 = [m1+2*sqrt(s21); flipdim(m1-2*sqrt(s21),1)];
% figure(3), fill([z; flipdim(z,1)], f1, [7 7 7]/8)
% hold on; plot(z, m1,'-r'); plot(x, y, '+k')

%% VFE - test 
disp('running vfe ...')
% init
hyp2.cov = hyp.cov;
hyp2.lik = hyp.lik;
hyp2.Xu = x(randperm(n,m),:);
covfunc1 = {@covSEiso};

D=size(x,2);

disp('training ...')
theta = vfeTrain([hyp.cov;hyp.lik;hyp2.Xu(:)],covfunc1,x,y,500);
hyp2.cov = theta(1:D+1);
hyp2.lik = theta(D+2);
hyp2.Xu = reshape(theta(D+3:end),m,D);

disp('predicting ...');
[mf, vf] = vfePredict(covfunc,hyp2,x,y,z);
vy = vf + exp(2*hyp2.lik);

fig1 = figure(10); cla, hold on;
PP = [0,0,16,10]; % *** paper position in centimeters
PS = PP(end-1:end); % paper size in centimeters
% Fonts
FontName = 'Times';
FSsm = 7; % small font size
FSmed = 10; % medium font size
FSlg = 11; % large font size
set(fig1, 'paperunits','centimeters','paperpositionmode','manual','paperposition', ...
    PP,'papersize',PS);
set(gcf, 'Color', 'None')
% So the figure is the same size on the screen as when it is printed:
pu = get(gcf,'PaperUnits');
pp = get(gcf,'PaperPosition');
set(gcf,'Units',pu,'Position',pp)
h2 = boundedline(z,mf,2*sqrt(vy),'r','alpha');
h1 = boundedline(z,m1,2*sqrt(s21),'b','alpha');
plot(hyp2.Xu,-1.8*ones(m,1),'ok')
h0 = plot(x, y, '+b');
xlabel('x','fontname',FontName,'FontSize',FSmed) % add axis labels
ylabel('y','fontname',FontName,'FontSize',FSmed)
legend([h1,h2],'Exact','VFE');
