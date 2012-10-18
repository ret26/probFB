function test_suite = test_addRamp
  initTestSuite;

% Tests: 
%
% Yramp = addRamp(Y,lens)
%

function test_alternative_computation

dispFigs=0;
T = 100;
Y = randn(T,1);
lens = 30;

Yramp1 = addRamp(Y,lens);

onSet = sin(2*pi*[1:lens]/(4*lens));
amp = [onSet,ones(1,T-2*lens),onSet(end:-1:1)]';

Yramp2 = Y.*amp;

if dispFigs ==1
  
  figure
  hold on
  plot(Y,'-k','linewidth',2)
  plot(Yramp1,'-r','linewidth',2)
  plot(Yramp1,'-b')
end

tol = 1e-6;
assertVectorsAlmostEqual(Yramp1,Yramp2,'absolute',tol,0)



function test_alternative_computation_2D

dispFigs=0;
T = 100;
Y = randn(T,2);
lens = [30;25];

Yramp1 = addRamp(Y,lens);
Yramp2 = zeros(size(Y));

for d=1:2
  onSet = sin(2*pi*[1:lens(d)]/(4*lens(d)));
  amp = [onSet,ones(1,T-2*lens(d)),onSet(end:-1:1)]';
  Yramp2(:,d) = Y(:,d).*amp;
end


tol = 1e-6;
assertVectorsAlmostEqual(Yramp1,Yramp2,'absolute',tol,0)


