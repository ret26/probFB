function test_suite = test_freq2probSpec
  initTestSuite;

% Tests: 
%
% function [om,lamx,varx] = freq2probSpec(fmax,df,varMa)
%

function testConsistent

dispFigs = 0;

fmax = [1/8,1/4]';
df = [1/32,1/80]';
mVar = [1,2]';

[om,lamx,varx] = freq2probSpec(fmax,df,mVar);
[fmax2,df2,mVar2] = probSpec2freq(om,lamx,varx);

tol = 1e-3;
assertVectorsAlmostEqual(fmax,fmax2,'absolute',tol,0)
assertVectorsAlmostEqual(df,df2,'absolute',tol,0)
assertVectorsAlmostEqual(mVar,mVar2,'absolute',tol,0)



function test_plot_spectra

dispFigs = 1;

fmax = [1/8]';
df = [1/32]';
mVar = [1]';

[om,lamx,varx] = freq2probSpec(fmax,df,mVar);

freqs = linspace(0,1/2,10000);
spec = get_pSTFT_spec(freqs,lamx,varx,om);

%mn_freq = sum(freqs.*spec/sum(spec));
%var_freq = (freqs-mn_freq).^2.*spec/sum(spec);
%[mn_freq,fmax]

if dispFigs==1
figure
hold on
plot(freqs,spec,'-k')
plot(fmax*[1,1],[0,max(spec)],'-r')
plot(fmax+[df,-df],max(spec)*[1,1]/2,'-b')

end

[val,loc]=max(spec);
fmaxTrue = freqs(loc);
ind = abs(diff(spec>val/2));
dfTrue = freqs(logical(ind));

tol = 1e-1;
assertVectorsAlmostEqual(fmax,fmaxTrue,'absolute',tol,0)
assertVectorsAlmostEqual(fmaxTrue-dfTrue(1),df,'absolute',tol,0)
assertVectorsAlmostEqual(dfTrue(2)-fmaxTrue,df,'absolute',tol,0)
