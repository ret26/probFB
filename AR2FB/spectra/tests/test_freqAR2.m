
function test_suite = test_freqAR2
  initTestSuite;

  %
  
function test_consistentWith_AR22freq

dispFigs = 0;

fmax = [1/8,1/4]';
df = [1/32,1/80]';
mVar = [1,2]';

[lamx,varx,dft] = freq2AR2(fmax,df,mVar,1);
[fmax2,df2, mVar2] = AR22freq(lamx,varx);

tol = 1e-5;
assertVectorsAlmostEqual(fmax,fmax2,'absolute',tol,0)
assertVectorsAlmostEqual(df,df2,'absolute',tol,0)
assertVectorsAlmostEqual(mVar,mVar2,'absolute',tol,0)

