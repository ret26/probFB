function test_denoiseSMGP_freq
%initTestSuite;
test_denoise_notraining_nonoise;
%test_denoise_notraining_withnoise;
end

function test_denoise_notraining_nonoise

varx = 0.5;
lenx = 100;
frex = 0.3;
T = 10000;
y = sampleGPSM(varx,lenx,frex,T);

% params = [varx lenx frex 0.00001]';

[mf,vf] = denoiseSMGP_freq(varx,lenx,frex,0.0001,y);

figure, hold on;
plot(1:T,y,'-b',1:T,mf,'-r');
legend('true','predicted');
xlabel('time'); ylabel('y')

%tol = 0.1;
%assertVectorsAlmostEqual(y(500:end-500),mf(500:end-500),'absolute',tol,0)
dl = mean((y-mf).^2)/var(y);
fprintf('smse %f\n',dl);

end

function test_denoise_notraining_withnoise

varx = 0.5;
lenx = 30;
frex = 0.3;
vary = 0.4;
T = 10000;
ytrue = sampleGPSM(varx,lenx,frex,T);
y = ytrue + 0.01*randn(T,1);

% params = [varx lenx frex vary]';

[mf,vf] = denoiseSMGP_freq(varx,lenx,frex,vary,y);

figure, hold on;
plot(1:T,y,'-b',1:T,ytrue,'-r',1:T,mf,'-g',1:T,mf+2*sqrt(vf),'--g',1:T,mf-2*sqrt(vf),'--g');
legend('noisy','true','predicted');
xlabel('time'); ylabel('y')

%tol = 0.1;
%assertVectorsAlmostEqual(ytrue(500:end-500),mf(500:end-500),'absolute',tol,0)
dl = mean((ytrue-mf).^2)/var(ytrue);
fprintf('smse %f\n',dl);
end