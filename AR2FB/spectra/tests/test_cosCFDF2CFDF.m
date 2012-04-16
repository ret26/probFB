
function test_suite = test_cosCFDF2CFDF
  initTestSuite;

  % tests both 
  % [CF,DF] = cosCFDF2CFDF(cosCF,cosDF)
  % and 
  % [cosCF,cosDF] = CFDF2cosCFDF(CF,DF)
  
function test_consistent

rand('state',6);
D = 5;
CF = rand(D,1)/2;
maxDF = min([1/2-CF,CF]')'*2;

DF = maxDF.*rand(D,1);

CF = 100/16000;
DF = 400/16000;

[cosCF,cosDF] = CFDF2cosCFDF(CF,DF);

% figure
% hold on
% f=linspace(0,1/2,100);
% plot(f,cos(2*pi*f),'-g')
% plot([CF(1),CF(1)],[0,cosCF(1)],'-k')
% plot(acos(cosCF(1)-cosDF(1)/2)*[1,1]/(2*pi),[0,cosCF(1)-cosDF(1)/2],'-r')
% plot(acos(cosCF(1)+cosDF(1)/2)*[1,1]/(2*pi),[0,cosCF(1)+cosDF(1)/2],'-b')

[CF2,DF2] = cosCFDF2CFDF(cosCF,cosDF);

% [CF,CF2]
% [DF,DF2]

%keyboard

tol = 1e-5;
assertVectorsAlmostEqual(CF,CF2,'absolute',tol,0)
assertVectorsAlmostEqual(DF,DF2,'absolute',tol,0)


