function test_suite = test_sampleGPSM
initTestSuite;
%test_sample1
end

function test_correctAutoCorrelation
dispFigs=1;
%randn('state',1)
varx = [1.5 0.4];
lenx = [5.32 8.9];
freqx = [0.02 0.3];

tau = 30;
T = 50000; % this needs to be large

y = sampleGPSM(varx,lenx,freqx,T);

K = length(varx);
autoCor = zeros(1,tau);
for k = 1:K
    autoCor = autoCor + varx(k)*exp(-1/(2*lenx(k)^2)*((0:tau-1)).^2)...
        .*cos(2*pi*freqx(k)*(0:tau-1));
end

autoCorEmp = zeros(tau,1);

for t=1:tau
    autoCorEmp(t) = mean(y(1:T-t+1).*y(t:T));
end

if dispFigs==1
    figure;
    hold on
    plot(autoCor,'-k','linewidth',2)
    plot(autoCorEmp,'-r','linewidth',1)
end

tol = 1e-1;
assertVectorsAlmostEqual(autoCor,autoCorEmp','absolute',tol,0)
end

function test_sample1
%{
    0.5000
  500.0000
    0.0400
    0.0001
    
   0.1331
   38.4573
    0.1004
    0.0002
%}

T = 10000;
y1 = sampleGPSM(0.5,500,0.04,T);
y2 = sampleGPSM(0.1331,38.4573,0.1004,T);
figure, plot(1:T,y1,1:T,y2)

end