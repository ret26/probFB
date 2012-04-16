
function test_suite = test_freq2AR2
  initTestSuite;

  % tests both freq2AR2 and AR22freq
  
function test_consistentWith_AR22freq

dispFigs = 0;

fmax = [1/8,1/4]';
df = [1/32,1/80]';
mVar = [1,2]';

[lamx,varx] = freq2AR2(fmax,df,mVar,1);
[fmax2,df2, mVar2] = AR22freq(lamx,varx);

tol = 1e-3;
assertVectorsAlmostEqual(fmax,fmax2,'absolute',tol,0)
assertVectorsAlmostEqual(df,df2,'absolute',tol,0)
assertVectorsAlmostEqual(mVar,mVar2,'absolute',tol,0)


function test_examples
% mainly for visualisation
dispFigs = 1;

fmax = 100/16000;
df = 140/16000;
mVar = 1;

[lamx,varx] = freq2AR2(fmax,df,mVar,1);

NumFreqs = 1000;
RngFreqs = [0,1/2];
[Freqs,Spec,fMAX,SpecMAX,dF1,dF2] = getSpecAR2(lamx,varx,NumFreqs, ...
						    RngFreqs);
if dispFigs==1
figure
hold on
plot(Freqs,Spec)
plot([fMAX,fMAX],[0,SpecMAX],'-r')
plot([dF1,dF2],[SpecMAX/2,SpecMAX/2],'-r')

disp(['desired DF: ',num2str(df),' attained DF: ',num2str(abs(dF1-dF2))])
end

tol =1e-3;
assertVectorsAlmostEqual(df,abs(dF1-dF2),'absolute',tol,0)