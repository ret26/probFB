% 28/7: 10 components
close all;
clear all;
clc
listfile = 'allfilelist.txt';
fid = fopen(listfile);
files = textscan(fid,'%s','delimiter','\n');
fclose(fid);
fdir = './datasets/audio/timit/';
myformat = 'fname: %s old_snr: %.2f new_snr: %.2f\n';
fname = [fdir files{1}{1} '.wav'];
[y,fs] = audioread(fname);
% normalise signal
rerate = 4;
y = resample(y,1,rerate);
y = y-mean(y);
y = y/std(y);
tau2 = 5;
tau1 = 100;
noEvals = 400;
T = length(y);
K = floor(T/tau1);
K = floor(K/4);
X = (1:tau1*K)'; % ignore samples at the end for now
Y = y(X);
noComps = 20;
window = 20;
%noComps = 20;
%window = 150;
overlap = 0.5;
params = initSECosParams(Y,fs/rerate,noComps,window,overlap,1);
drawnow;
keyboard
wnoise = randn(tau1*K,1);
noiselevel = 0;
Ynoisy = addnoise(Y,wnoise,noiselevel);
missing = zeros(size(X)) == 1;
missingStack = reshape(missing,[tau1,K])';
covfunc = {@covSEiso};
theta_init = log(params);
[theta_end,nlml] = testing_trainSECosMix(theta_init,noComps,covfunc,...
    X,Ynoisy,tau1,tau2,missingStack,noEvals);
[fest,vest] = testing_predictSECosMix(theta_end,noComps,covfunc,...
    X,Ynoisy,tau1,tau2,missingStack);
old_snr = snr(Y,Y-Ynoisy);
new_snr = snr(Y,Y-fest);
fprintf(myformat,fname,old_snr,new_snr);
