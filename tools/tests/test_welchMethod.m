
function test_suite = test_welchMethod
  initTestSuite;


function test_sinusoid

dispFigs=0;
T = 10000;

Ts = [30,50,100];
%as = [1/2,1/4,3/4];
as = [1,1,1];

y = as(1)*sin(2*pi*[1:T]'/Ts(1));

K = length(Ts);
for k=2:K
  y = y+as(k)*sin(2*pi*[1:T]'/Ts(k));
end

numFreq = 400;
ovLp = 10;
[pg,varpg] = welchMethod(y,numFreq,ovLp);

freq = linspace(0,1/2,numFreq);

  cumSumSpecTr = zeros(numFreq,1);
  for k=1:K
    cumSumSpecTr = cumSumSpecTr+as(k).^2*(freq'>1./Ts(k))/4;
  end

if dispFigs==1
  figure
  subplot(2,1,1)
  hold on
  for k=1:K
    plot(1./Ts(k)*[1,1],[0,as(k).^2]/4,'-r','linewidth',2)
  end
  plot(freq,pg,'-k')
  
  subplot(2,1,2)
  hold on
  plot(freq,cumSumSpecTr,'-r','linewidth',2)
  plot(freq,cumsum(pg),'-k')
end

tol = 3e-1;
assertVectorsAlmostEqual(cumsum(pg),cumSumSpecTr,'absolute',tol,0)

function test_AR

dispFigs=1;
T = 100000;

fmax = 1/8;
df = 1/32;
varMa = 1;

[lamx,varx,dft] = freq2AR2(fmax,df,varMa);

[y,x] = sampleAR2FB(lamx,varx,0,T);

numFreq = 200;
ovLp = 195;
[pg,varpg] = welchMethod(y,numFreq,ovLp);

freq = linspace(0,1/2,numFreq);

[Freqs,specTr,fMAX,SpecMAX,dF1,dF2] = getSpecAR2(lamx,varx,numFreq,[0,1/2]);

specTr = 1/2*specTr/sum(specTr);

cumSumSpecTr = cumsum(specTr);


if dispFigs==1
  figure
  subplot(2,1,1)
  hold on
  tol = 3;
  up = pg+tol*sqrt(varpg);
  low = pg-tol*sqrt(varpg);
  EBx = [freq(:);freq(end:-1:1)'];
  EBy = [up;low(end:-1:1)];
  
  patch(EBx,EBy,[1,0.5,0.5],'edgecolor','none');
  plot(freq,pg,'-k')
  plot(Freqs,specTr,'-b')

  subplot(2,1,2)
  hold on
  plot(freq,cumSumSpecTr,'-b','linewidth',2)
  plot(freq,cumsum(pg),'-k')
end

tol = 3e-1;
assertVectorsAlmostEqual(cumsum(pg),cumSumSpecTr','absolute',tol,0)