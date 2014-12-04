
close all;
clear all;
clc
listfile = '../../../tsgp/datasets/audio/timit/allfilelist.txt';
fid = fopen(listfile);
files = textscan(fid,'%s','delimiter','\n');
fclose(fid);
fdir = '../../../tsgp/datasets/audio/timit/';
myformat = 'fname: %s old_snr: %.2f new_snr: %.2f\n';
fname = [fdir files{1}{1} '.wav'];
[y,fs] = audioread(fname);
% normalise signal
rerate = 4;
y = resample(y,1,rerate);
y = y-mean(y);
y = y/std(y);
y = y(800:2000);
tau2 = 3;
tau1 = 50;
noEvals = 400;
T = length(y);
K = floor(T/tau1);
M = 5;
wnoise = randn(T,1);
noiselevel = 40;
ynoisy = addnoise(y,wnoise,noiselevel);
% keyboard
noComps = 10;
window = 3;
params_learnt = [];
m_sig_est = [];
v_sig_est = [];
for k = 1:K
    fprintf('k=%i/%i\n',k,K)
    kstart = k-M;
    if kstart<0; kstart = 0; end
    indWindow = (kstart*tau1+1):k*tau1;
    indFrame  = (k-1)*tau1 + (1:tau1);
    X = indWindow';
    %Y = y(X);
    Ynoisy = ynoisy(X);
    overlap = 0.5;
    params = initSECosParams(Ynoisy,fs/rerate,noComps,window,overlap,0);
    %drawnow;
    missing = zeros(size(X)) == 1;
    missingStack = reshape(missing,[tau1,k-kstart])';
    covfunc = {@covSEiso};
    weights = 0.8.^(k-kstart-1:-1:0);
    theta_init = log(params);
    [theta_end,nlml] = testing_trainAdaptive(theta_init,noComps,covfunc,...
        X,Ynoisy,tau1,tau2,missingStack,noEvals,weights);
    [fest,vest] = testing_predictAdaptive(theta_end,noComps,covfunc,...
        X,Ynoisy,tau1,tau2,missingStack);
    params_learnt = [params_learnt theta_end];
    m_sig_est = [m_sig_est;fest];
    v_sig_est = [v_sig_est;vest];
end
Y = y(1:tau1*K);
Ynoisy = ynoisy(1:tau1*K);
old_snr = snr(Y,Y-Ynoisy);
new_snr = snr(Y,Y-m_sig_est);
fprintf(myformat,fname,old_snr,new_snr);

%%
figure(5);
subplot(3,1,1), plot(Y)
subplot(3,1,2), plot(Ynoisy)
subplot(3,1,3), plot(m_sig_est)
%%
figure(6);
subplot(3,1,1), plot(Y)
subplot(3,1,2), plot(Ynoisy-Y)
subplot(3,1,3), plot(m_sig_est-Y)
%%
freqs = exp(params_learnt(2*noComps+1:end-1,:));
figure(7), plot(freqs','ob'); ylim([0,0.5]);
